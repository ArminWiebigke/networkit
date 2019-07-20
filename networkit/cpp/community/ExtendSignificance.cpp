/*
 * ExtendSignificance.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#include <set>

#include "ExtendSignificance.h"
#include "../auxiliary/Timer.h"
#include "../auxiliary/ParseString.h"
#include "../coarsening/ParallelPartitionCoarsening.h"
#include "../oslom/Stochastics.h"
#include "../structures/AdjacencyArray.h"
#include "EgoSplitting.h"

#define W(x) #x << "=" << x << ", "
#define true_or_throw(cond, msg) if (!cond) throw std::runtime_error(msg)

namespace NetworKit {

ExtendSignificance::ExtendSignificance(const EgoNetData &egoNetData, const Partition &basePartition,
                                       count maxCandidates)
		: ExtendScore(egoNetData, maxCandidates),
		  basePartition(basePartition),
		  sigTable(egoNetData.sigTable),
		  useSigMemo(parameters.at("useSigMemo") == "Yes"),
		  mergeGroups(parameters.at("signMerge") == "Yes"),
		  sortGroupsStrat(parameters.at("sortGroups") == "significance"),
		  maxSignificance(Aux::stringToDouble(parameters.at("maxSignificance"))),
		  maxGroupCnt(std::stoi(parameters.at("maxGroupsConsider"))),
		  minEdgesToGroup(std::stoi(parameters.at("minEdgesToGroupSig"))) {
	double maxSig = maxSignificance;
	// For a given number of external nodes, calculate the minimal r-score that is significant
	// We can then simply compare the r-values instead of calculating the ordered statistics
	// Only works if the position is known before. 0 ( = best random node) in this case.
	sigTable.trySetValueFunc([maxSig](index extNodes) {
		double minRScore = 0.125;
		double step = 0.0625;
		int iterations = 20;
		for (int i = 0; i < iterations; ++i) {
			double sig = Stochastics::order_statistics_left_cumulative(extNodes, extNodes,
			                                                           minRScore);
			if (sig < maxSig) {
				minRScore += step;
			} else {
				minRScore -= step;
			}
			step /= 2;
		}
		return minRScore;
	});
	significantGroup.resize(G.upperNodeIdBound(), none); // TODO: only do this once
}

void ExtendSignificance::run() {
	Aux::Timer timer;
	timer.start();
	Stochastics::init(G.numberOfEdges());

	// Get coarse graph from base partition
	createCoarseGraph();
	addTime(timer, "0    Coarsening");

	// Get candidates
	edgesToGroups.resize(G.upperNodeIdBound()); // TODO: only do this once
	addTime(timer, "2    Init candidates");
	std::vector<node> candidates = getCandidates();
	addTime(timer, "3    Find candidates");

	// Insert a node that represents all external nodes and add the outgoing edges
	addExtEdges(candidates);
	addTime(timer, "4    Calc outgoing edges");
	groupStubs = calcGroupStubsCounts();
	addTime(timer, "4a    Calculate stub counts");

	// Sort candidates by number of edges into the egoNet, filter bad candidates (less than 3 edges)
	candidatesSorted = sortCandidatesByEdges(candidates); // TODO: Why do we have to sort?
	addTime(timer, "5    Sort candidates");

	if (parameters.at("onlyCheckSignOfMaxCandidates") == "Yes") {
		double evalFactor = Aux::stringToDouble(parameters.at("evalSignFactor"));
		count evalCandidates = evalFactor * maxExtendedNodes;
		if (candidatesSorted.size() > evalCandidates)
			candidatesSorted.resize(evalCandidates);
	}


//	for (node v : directNeighbors)
//			directedG.forEdgesOf(v, countEdge);
//	for (node w : candidates) {
//		directedG.forEdgesOf(w, [&](node, node v, edgeweight weight) {
//			if (isDirectNeighbor(v))
//				edgesToGroups[w][egoToCoarse[egoMapping.local(v)]] += weight;
//		});
//	}
//	auto candidatesDir = candidates;
//	candidates.clear();
//	auto edgeScoresDir = edgesToGroups;
//	edgesToGroups = std::vector<std::map<node, double>>(G.upperNodeIdBound());
//	for (node v : directNeighbors) {
//		G.forEdgesOf(v, countEdge);
//	}
//	assert(candidatesDir.size() <= candidates.size());
//	for (node w : candidates){
//		for (auto pair : edgesToGroups[w])
//			assert(edgeScoresDir[w][pair.first] <= edgesToGroups[w][pair.first]);
//	}

	// Check
	auto check = [&]() {
		for (node v : candidates) {
			std::unordered_map<node, count> groupEdges;
			G.forNeighborsOf(v, [&](node w) {
				if (egoMapping.isMapped(v))
					groupEdges[egoToCoarse[w]] += 1;
			});
			for (auto pair : groupEdges)
				assert(pair.second == coarseGraph.weight(v, pair.first));
		}
		return true;
	};
	assert(check());

	// Calculate significances
	checkCandidates("6");
	addTime(timer, "6    Calc significance one");

	// Calculate significances, allowing slightly worse candidates (based on the number of already
	// found significant candidates)
	// TODO: Just store the best significance, then test against new ordered statistic.
	//  Need to store significance and number of external nodes
	int iterations = std::stoi(parameters.at("secondarySigExtRounds"));
	for (int i = 0; i < iterations; ++i) {
		secondRound();
	}
	addTime(timer, "8    Calc significance two");

#ifndef NDEBUG
	std::set<node> resultSet;
	for (auto pair : result) {
		node v = pair.first;
		assert(resultSet.count(v) == 0);
		resultSet.insert(v);
	}
	assert(resultSet.size() == result.size());
#endif

	hasRun = true;
}

void ExtendSignificance::updateCandidates() {
	for (node w : addedCandidates) {
		node group = significantGroup[w];
		// Update group stubs
		// TODO: edges between added candidates
		groupStubs.groupTotal[group] += G.degree(w);
		groupStubs.groupOutgoing[group] += G.degree(w) - 2 * edgesToGroups[w][group];
		groupStubs.externalNodes[group] -= 1;
		groupStubs.externalStubs[group] -= G.degree(w);
		G.forNeighborsOf(w, [&](node v) {
			// Update other candidates
			if (!edgesToGroups[v].empty()) {
				// v is a neighbor of neighbor
				edgesToGroups[v][group] += 1;
			}
			if (significantGroup[v] == group) {
				// v is in the same group
				count newInternalStubs = 2;
				if (addedCandidates.count(v)) {
					// Don't count this edge twice if we add both candidats at the same time
					newInternalStubs = 1;
				}
				groupStubs.groupOutgoing[group] -= newInternalStubs;
			}
		});
		edgesToGroups[w].clear();
	}
}

void ExtendSignificance::createCoarseGraph() {
	ParallelPartitionCoarsening coarsening(egoGraph, basePartition);
	coarsening.run();
	coarseGraph = coarsening.getCoarseGraph();
	coarseToEgo = coarsening.getCoarseToFineNodeMapping();
	egoToCoarse = coarsening.getFineToCoarseNodeMapping();
	coarseSizes.resize(coarseGraph.upperNodeIdBound());
	coarseGraph.forNodes([&](node p) {
		coarseSizes[p] = coarseToEgo[p].size();
	});
	numGroups = coarseGraph.upperNodeIdBound();
}

void ExtendSignificance::secondRound() {
	scorePenalty = 0.5;
	std::string strat = parameters.at("sigSecondRoundStrat");
	if (strat == "updateCandidates") {
		updateCandidates();
	} else if (strat == "orderStat") {
		orderStatPos = std::floor(Aux::stringToDouble(parameters.at("orderedStatPos")) * result.size());
		if (orderStatPos == 0)
			return;
	} else {
		throw std::runtime_error(strat + " is not a valid strategy for the second signif. round!");
	}
	removeAddedAsCandidates();
	checkCandidates("8");
}

void ExtendSignificance::removeAddedAsCandidates() {
	auto newEnd = std::remove_if(candidatesSorted.begin(), candidatesSorted.end(),
	                             [&](std::pair<count, node> scorePair) {
		                             return addedCandidates.count(scorePair.second) > 0;
	                             });
	candidatesSorted.resize(newEnd - candidatesSorted.begin());
	addedCandidates.clear();
}

std::vector<std::pair<count, node>>
ExtendSignificance::sortCandidatesByEdges(const std::vector<node> &candidates) const {
	std::vector<std::pair<count, node>> sorted;
	sorted.reserve(candidates.size());
	count minEdges = 3;
	for (node w : candidates) {
		count edgesIntoEgo = std::accumulate(edgesToGroups[w].begin(), edgesToGroups[w].end(),
		                                     (count) 0);
		if (edgesIntoEgo >= minEdges)
			sorted.emplace_back(edgesIntoEgo, w);
	}
	std::sort(sorted.rbegin(), sorted.rend());
	return sorted;
}

node
ExtendSignificance::addExtEdges(std::vector<node> &candidates) {
	node externalNode = coarseGraph.addNode();
	assert(numGroups == externalNode);
	coarseSizes.push_back(G.numberOfNodes() - egoGraph.numberOfNodes());
	// Calc outgoing edges for each group
	candidates.push_back(egoNode); // Add edges to egoNode, but remove it as candidate later
	for (node w : candidates) {
		for (node p = 0; p < edgesToGroups[w].size(); ++p) {
			count numEdges = edgesToGroups[w][p];
			coarseGraph.increaseWeight(p, externalNode, numEdges);
		}
	}
	candidates.pop_back();
	return externalNode;
}

std::vector<node>
ExtendSignificance::getCandidates() {
	std::vector<node> candidates;
	std::vector<node> directNeighbors = egoMapping.globalNodes();
	auto isDirectNeighbor = [&](node x) {
		return egoMapping.isMapped(x);
	};
	auto countEdge = [&](node vLoc, node w, edgeweight weight) {
		if (edgesToGroups[w].empty()) {
			candidates.push_back(w);
			edgesToGroups[w].resize(numGroups);
		}
		edgesToGroups[w][vLoc] += (count) weight;
	};
	auto removeNeighborCandidates = [&]() {
		for (node w : candidates) {
			if (isDirectNeighbor(w)) {
				edgesToGroups[w].clear();
			}
		}
		auto newEnd = std::remove_if(candidates.begin(), candidates.end(), [&](node w) {
			return (isDirectNeighbor(w));
		});
		candidates.resize(newEnd - candidates.begin());
	};

	if (parameters.at("extendOverDirected") == "Yes") {
		for (node v : directNeighbors) {
			node vLoc = egoToCoarse[egoMapping.local(v)];
			directedG.forEdgesOf(v, [&](node, node w, edgeweight weight) {
				countEdge(vLoc, w, weight);
			});
		}
		removeNeighborCandidates();
		if (parameters.at("extendDirectedBack") == "Yes") {
			for (node w : candidates) {
				directedG.forEdgesOf(w, [&](node, node v, edgeweight weight) {
					if (isDirectNeighbor(v))
						edgesToGroups[w][egoToCoarse[egoMapping.local(v)]] += (count) weight;
				});
			}
		}
	} else {
		for (node v : directNeighbors) {
			node vLoc = egoToCoarse[egoMapping.local(v)];
			G.forEdgesOf(v, [&](node, node w, edgeweight weight) {
				countEdge(vLoc, w, weight);
			});
		}
		removeNeighborCandidates();
	}
	removeEgoNodeCandidate(candidates);
	return candidates;
}

void ExtendSignificance::removeEgoNodeCandidate(std::vector<node> &candidates) {
	for (int i = 0; i < candidates.size(); ++i) {
		if (candidates[i] == egoNode) {
			candidates[i] = candidates.back();
			candidates.pop_back();
			break;
		}
	}
	edgesToGroups[egoNode].clear();
}

void
ExtendSignificance::checkCandidates(const std::string &t_prefix) {
	Aux::Timer timer;
	timer.start();

	for (const auto &pair : candidatesSorted) {
		if (enoughCandidates())
			break;
		node v = pair.second;
		checkCandidate(t_prefix, timer, v);
	}
}

bool ExtendSignificance::enoughCandidates() const { return result.size() >= maxExtendedNodes; }

void
ExtendSignificance::checkCandidate(const std::string &t_prefix, Aux::Timer &timer, node v) {
	std::vector<std::pair<double, node>> groupEdges = sortGroupsByEdges(v);
	if (groupEdges.empty())
		return;
	addTime(timer, t_prefix + "b    Sort groups by edge count");

	std::vector<double> groupSigs(numGroups);
	bool added = checkSingleGroups(v, groupEdges, groupSigs);
	addTime(timer, t_prefix + "c    Significance single");
	if (added)
		return;

	if (!mergeGroups)
		return;
	checkMergedGroups(t_prefix, timer, v, groupEdges, groupSigs);
	addTime(timer, t_prefix + "f    Significance merge");
}

std::vector<std::pair<double, node>> ExtendSignificance::sortGroupsByEdges(node v) const {
	std::vector<std::pair<double, node>> groupEdges;
	groupEdges.reserve(edgesToGroups[v].size());
	for (node p = 0; p < edgesToGroups[v].size(); ++p) {
		count numEdges = edgesToGroups[v][p];
		if (numEdges >= minEdgesToGroup)
			groupEdges.emplace_back(numEdges, p);
	}
	std::sort(groupEdges.rbegin(), groupEdges.rend());
	return groupEdges;
}

GroupStubs ExtendSignificance::calcGroupStubsCounts() const {
	GroupStubs stubs(numGroups);
	for (count p = 0; p < numGroups; ++p) {
		if (!coarseGraph.hasNode(p))
			continue;
		auto groupDegree = (count) coarseGraph.weightedDegree(p);
		stubs.groupTotal[p] = groupDegree;
		stubs.groupOutgoing[p] = groupDegree - 2 * (count) coarseGraph.weight(p, p);
		stubs.externalStubs[p] = 2 * G.numberOfEdges() - groupDegree;
		stubs.externalNodes[p] = G.numberOfNodes() - coarseSizes[p];
	}
	return stubs;
}

bool ExtendSignificance::checkSingleGroups(node v,
                                           const std::vector<std::pair<double, node>> &groupEdges,
                                           std::vector<double> &groupSigs) {
	count nodeDegree = G.degree(v);
	bool added = false;
	count calcedGroups = 0;
	for (auto it : groupEdges) {
		if (calcedGroups >= maxGroupCnt)
			break;
		++calcedGroups;
		count numEdges = (int) it.first;
		node p = it.second;
//			addTime(timer, t_prefix + "a    Loop b");

		double significance = calcScore(nodeDegree, (int) numEdges, groupStubs.groupOutgoing[p],
		                                groupStubs.externalStubs[p], groupStubs.externalNodes[p]);
		groupSigs[p] = significance;
		added = addIfSignificant(v, significance, p);
//			addTime(timer, t_prefix + "e    Add if significant");
		if (added)
			break;
	}
	return added;
}

bool ExtendSignificance::addIfSignificant(node v, double significance, node group) {
	if (significance <= maxSignificance) {
		double score = 1 - significance - scorePenalty;
		result.emplace_back(v, score);
		significantGroup[v] = group;
		addedCandidates.insert(v);
		return true;
	}
	return false;
};

bool
ExtendSignificance::checkMergedGroups(const std::string &t_prefix, Aux::Timer &timer, node v,
                                      std::vector<std::pair<double, node>> &groupEdges,
                                      const std::vector<double> &groupSigs) {
	count nodeDegree = G.degree(v);
	if (sortGroupsStrat) {
		std::sort(groupEdges.begin(), groupEdges.end(), [&](std::pair<double, node> a,
		                                                    std::pair<double, node> b) {
			return (groupSigs[a.second]) > groupSigs[b.second];
		});
	}
	addTime(timer, t_prefix + "d    Sort groups by significance");

	auto it = groupEdges.begin();
	node bestGroup = it->second;
	std::vector<node> mergedGroups{bestGroup};
	count totalStubsMerged = groupStubs.groupTotal[bestGroup];
	count outStubsMerged = groupStubs.groupOutgoing[bestGroup];
	count extStubsMerged = groupStubs.externalStubs[bestGroup];
	count extNodesMerged = groupStubs.externalNodes[bestGroup];
	count numEdgesMerged = (int) it->first;
	bool added = false;

	for (++it; it < groupEdges.end(); ++it) {
		node group = it->second;
		if (groupSigs[group] == 0.0)
			break;
		count newInternalStubs = 0;
		for (node w : mergedGroups) {
			newInternalStubs += 2 * (count) coarseGraph.weight(w, group);
		}
		// Groups have to be connected (?)
		if (newInternalStubs == 0)
			continue;
		// Merge group
		extStubsMerged -= groupStubs.groupTotal[group];
		totalStubsMerged += groupStubs.groupTotal[group];
		outStubsMerged += groupStubs.groupOutgoing[group] - newInternalStubs;
		extNodesMerged -= coarseSizes[group];
		numEdgesMerged += (count) it->first;
		mergedGroups.push_back(group);
//			addTime(timer, t_prefix + "g    Merge groups");

		// Calculate new significance
		double significance = calcScore(nodeDegree, numEdgesMerged, outStubsMerged,
		                                extStubsMerged, extNodesMerged);
		added = addIfSignificant(v, significance, bestGroup);
//			addTime(timer, t_prefix + "h    Add if significant (merge)");
		if (added) {
//				std::cout << mergedGroups.size() << " merged groups " << nodeScores.back().second
//				          << " (" << numEdgesMerged << "/" << nodeDegree << ") -> "
//				          << G.numberOfNodes() - extNodesMerged
//				          << std::endl;
			break;
		}
	}
	return added;
}

double
ExtendSignificance::calcScore(int nodeDegree, int kIn, int grOut, int groupExtStubs,
                              int extNodes) const {
	// TODO: bootInterval anschauen
	count openStubs = groupExtStubs + grOut -
	                  nodeDegree; // Remove node stubs from the free external stubs (?)
//	std::cout << W(grOut) << W(groupExtStubs) << W(nodeDegree) << W(openStubs) << std::endl;
	double rScore = Stochastics::compute_simple_fitness(kIn, grOut, openStubs, nodeDegree);

	// TODO: Possible speedup: for each node degree (and group), store the largest kIn that
	//   was not significant. If the current kIn is smaller than the stored one, we know that
	//   the current candidate is not significant
	//   Could do the same in the other direction, store the smalles kIn that was significant

	double returnVal = 1.0;
	if (useSigMemo && orderStatPos == 0) {
		// TODO: Why does this add more nodes?
		double minR = sigTable.getValue(extNodes);
		if (rScore <= minR * 0.99) {
			returnVal = rScore;
			assert(Stochastics::order_statistics_left_cumulative(
					extNodes, extNodes - orderStatPos, rScore) <= maxSignificance);
		}
		if (rScore >= minR * 1.01) {
			assert(Stochastics::order_statistics_left_cumulative(
					extNodes, extNodes - orderStatPos, rScore) >= maxSignificance);
		}

	} else {
		returnVal = Stochastics::order_statistics_left_cumulative(
				extNodes, extNodes - orderStatPos, rScore);

	}
	return returnVal;
}

std::string ExtendSignificance::toString() const {
	return "ExtendSignificance";
}

bool ExtendSignificance::isParallel() const {
	return false;
}

} /* namespace NetworKit */