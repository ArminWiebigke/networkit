/*
 * ExtendSignificance.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#include <set>

#include "ExtendSignificance.h"
#include "../auxiliary/Timer.h"
#include "../coarsening/ParallelPartitionCoarsening.h"
#include "../oslom/Stochastics.h"
#include "../structures/AdjacencyArray.h"
#include "EgoSplitting.h"

namespace NetworKit {

ExtendSignificance::ExtendSignificance(const EgoNetData &egoNetData,
                                       const Partition &basePartition)
		: ExtendScore(egoNetData), basePartition(basePartition), sigTable(egoNetData.sigTable),
		  subtractNodeDegree(parameters.at("subtractNodeDegree") == "Yes"),
		  useSigMemo(parameters.at("useSigMemo") == "Yes"),
		  mergeGroups(parameters.at("signMerge") == "Yes"),
		  sortGroupsStrat(parameters.at("sortGroups") == "significance"),
		  maxSignificance(std::stod(parameters.at("maxSignificance"))),
		  maxGroupCnt(std::stoi(parameters.at("maxGroupsConsider"))),
		  minEdgesToGroup(std::stoi(parameters.at("minEdgesToGroupSig"))) {
	double maxSig = maxSignificance;
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
		std::cout << extNodes << ": " << minRScore
		          << std::endl;
		return minRScore;
	});
}

void ExtendSignificance::run() {
	Aux::Timer timer;
	timer.start();
	Stochastics::init(G.numberOfEdges());

	// Get coarse graph from base partition
	createCoarseGraph();
	//	addTime(timer, "0    Coarsening");

	// Get candidates
	edgeScores.resize(G.upperNodeIdBound());
	addTime(timer, "2    Init candidates");
	std::vector<node> candidates = getCandidates();
	addTime(timer, "3    Find candidates");

	// Insert a node that represents all external nodes and add the outgoing edges
	node externalNode = addExtEdges(candidates);
	addTime(timer, "4    Calc outgoing edges");

	// Sort candidates by number of edges into the egoNet
	std::vector<std::pair<count, node>> candidatesSorted = sortCandidates(candidates);
	addTime(timer, "5    Sort candidates");

//	for (node v : directNeighbors)
//			directedG.forEdgesOf(v, countEdge);
//	for (node w : candidates) {
//		directedG.forEdgesOf(w, [&](node, node v, edgeweight weight) {
//			if (isDirectNeighbor(v))
//				edgeScores[w][egoToCoarse[egoMapping.local(v)]] += weight;
//		});
//	}
//	auto candidatesDir = candidates;
//	candidates.clear();
//	auto edgeScoresDir = edgeScores;
//	edgeScores = std::vector<std::map<node, double>>(G.upperNodeIdBound());
//	for (node v : directNeighbors) {
//		G.forEdgesOf(v, countEdge);
//	}
//	assert(candidatesDir.size() <= candidates.size());
//	for (node w : candidates){
//		for (auto pair : edgeScores[w])
//			assert(edgeScoresDir[w][pair.first] <= edgeScores[w][pair.first]);
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
	result = calcSignficance(externalNode, 0, "6", candidatesSorted);
	addTime(timer, "6    Calc significance one");

	// Calculate significances, allowing slightly worse candidates (based on the number of already
	// found significant candidates)
	// TODO: Just store the best significance, then test against new ordered statistic
	//  need to store significance and number of external nodes
	double orderedStatPosition = std::stod(parameters.at("orderedStatPos")) * result.size();
	if (orderedStatPosition > 0.0) {
		secondSigRound(externalNode, candidatesSorted, orderedStatPosition);
		addTime(timer, "8    Calc significance two");
	}

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
}

void ExtendSignificance::secondSigRound(node externalNode,
                                        std::vector<std::pair<count, node>> &candidatesSorted,
                                        double orderedStatPosition) {
	std::unordered_set<node> sigCandidates;
	for (auto pair : result) {
		node v = pair.first;
		sigCandidates.insert(v);
	}
	auto newEnd = std::remove_if(candidatesSorted.begin(), candidatesSorted.end(),
	                             [&](std::pair<count, node> scorePair) {
		                             return sigCandidates.count(scorePair.second) > 0;
	                             });
	candidatesSorted.resize(newEnd - candidatesSorted.begin());

	auto scores2 = calcSignficance(externalNode,
	                               orderedStatPosition, "8", candidatesSorted);
	for (auto pair : scores2) {
		// TODO: make sure that this score is low and > 0
		result.emplace_back(pair.first, pair.second - 0.5);
	}
}

std::vector<std::pair<count, node>>
ExtendSignificance::sortCandidates(const std::vector<node> &candidates) const {
	std::vector<std::pair<count, node>> candidatesSorted;
	candidatesSorted.reserve(candidates.size());
	count minEdges = 3;
	for (node w : candidates) {
		count edgesIntoEgo = (int) std::accumulate(edgeScores[w].begin(), edgeScores[w].end(), 0);
		if (edgesIntoEgo >= minEdges)
			candidatesSorted.emplace_back(edgesIntoEgo, w);
	}
	std::sort(candidatesSorted.rbegin(), candidatesSorted.rend());
	return candidatesSorted;
}

node
ExtendSignificance::addExtEdges(std::vector<node> &candidates) {
	node externalNode = coarseGraph.addNode();
	coarseSizes.push_back(G.numberOfNodes() - egoGraph.numberOfNodes());
	// Calc outgoing edges for each group
	candidates.push_back(egoNode); // Add edges to egoNode, but remove it as candidate later
	for (node w : candidates) {
		for (node p = 0; p < edgeScores[w].size(); ++p) {
			double numEdges = edgeScores[w][p];
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
	count numGroups = coarseGraph.upperNodeIdBound();
	auto isDirectNeighbor = [&](node x) {
		return egoMapping.isMapped(x);
	};
	auto countEdge = [&](node vLoc, node w, edgeweight weight) {
		if (edgeScores[w].empty()) {
			candidates.push_back(w);
			edgeScores[w].resize(numGroups);
		}
		edgeScores[w][vLoc] += weight;
	};
	auto removeNeighborCandidates = [&]() {
		auto newEnd = std::remove_if(candidates.begin(), candidates.end(), [&](node w) {
			return (isDirectNeighbor(w));
		});
		candidates.resize(newEnd - candidates.begin());
//		addTime(timer, "4    Remove neighbors as candidates");
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
						edgeScores[w][egoToCoarse[egoMapping.local(v)]] += weight;
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

void ExtendSignificance::removeEgoNodeCandidate(std::vector<node> &candidates) const {
	for (int i = 0; i < candidates.size(); ++i) {
		if (candidates[i] == egoNode) {
			candidates[i] = candidates.back();
			candidates.pop_back();
			break;
		}
	}
}

std::vector<std::pair<node, double>>
ExtendSignificance::calcSignficance(
		node numGroups, double orderedStatPosition,
		const std::string &t_prefix,
		const std::vector<std::pair<count, node>> &candidatesSorted) const {
	Aux::Timer timer;
	timer.start();

	GroupStubs groupStubs = calcGroupStubsCounts(numGroups);
	addTime(timer, t_prefix + "5    Calculate stub counts");

	index statPosInt = (int) std::floor(orderedStatPosition);
	std::vector<std::pair<node, double>> nodeScores;
	for (const auto &pair : candidatesSorted) {
		node v = pair.second;
		checkForSignificance(numGroups, t_prefix, timer, groupStubs, statPosInt, nodeScores, v);
	}
	return nodeScores;
}

void
ExtendSignificance::checkForSignificance(node numGroups, const std::string &t_prefix,
                                         Aux::Timer &timer,
                                         const GroupStubs &groupStubs, index statPosInt,
                                         std::vector<std::pair<node, double>> &nodeScores,
                                         node v) const {
	std::vector<std::pair<double, node>> groupEdges = sortGroupsByEdges(v);
	if (groupEdges.empty())
		return;
	addTime(timer, t_prefix + "b    Sort groups by edge count");

	std::vector<double> groupSigs(numGroups);
	count calcedGroups = 0;
	bool added = checkSingleGroups(groupStubs, statPosInt, nodeScores, v, groupEdges, groupSigs,
	                               calcedGroups);
	addTime(timer, t_prefix + "c    Significance single");
	if (added)
		return;

	if (!mergeGroups)
		return;
	checkMergedGroups(t_prefix, timer, groupStubs, statPosInt, v, groupEdges,
	                  groupSigs, calcedGroups, nodeScores);
	addTime(timer, t_prefix + "f    Significance merge");
}

std::vector<std::pair<double, node>> ExtendSignificance::sortGroupsByEdges(node v) const {
	std::vector<std::pair<double, node>> groupEdges;
	groupEdges.reserve(edgeScores[v].size());
	for (node p = 0; p < edgeScores[v].size(); ++p) {
		double numEdges = edgeScores[v][p];
		if (numEdges >= minEdgesToGroup)
			groupEdges.emplace_back(numEdges, p);
	}
	std::sort(groupEdges.rbegin(), groupEdges.rend());
	return groupEdges;
}

GroupStubs ExtendSignificance::calcGroupStubsCounts(node numGroups) const {
	GroupStubs groupStubs(numGroups);
	for (count p = 0; p < numGroups; ++p) {
		if (!coarseGraph.hasNode(p))
			continue;
		count groupDegree = (int) coarseGraph.weightedDegree(p);
		groupStubs.groupTotal[p] = groupDegree;
		groupStubs.groupOutgoing[p] = groupDegree - 2 * (int) coarseGraph.weight(p, p);
		groupStubs.externalStubs[p] = G.numberOfEdges() * 2 - groupStubs.groupTotal[p];
		groupStubs.externalNodes[p] = G.numberOfNodes() - coarseSizes[p];
	}
	return groupStubs;
}

bool ExtendSignificance::checkSingleGroups(const GroupStubs &groupStubs, index statPosInt,
                                           std::vector<std::pair<node, double>> &nodeScores, node v,
                                           const std::vector<std::pair<double, node>> &groupEdges,
                                           std::vector<double> &groupSigs,
                                           count &calcedGroups) const {
	count nodeDegree = G.degree(v);
	nodeScores.reserve(groupEdges.size());
	bool added = false;
	for (auto it : groupEdges) {
		if (calcedGroups++ >= maxGroupCnt)
			break;
		count numEdges = (int) it.first;
		node p = it.second;
//			addTime(timer, t_prefix + "a    Loop b");

		double significance = calcScore(nodeDegree, (int) numEdges, groupStubs.groupOutgoing[p],
		                                groupStubs.externalStubs[p], groupStubs.externalNodes[p],
		                                statPosInt);
		groupSigs[p] = significance;
		added = addIfSignificant(nodeScores, v, significance);
//			addTime(timer, t_prefix + "e    Add if significant");
		if (added)
			break;
	}
	return added;
}

bool ExtendSignificance::addIfSignificant(std::vector<std::pair<node, double>> &nodeScores,
                                          node v, double significance) const {
	if (significance <= maxSignificance) {
		double score = 1 - significance;
		nodeScores.emplace_back(v, score);
		return true;
	}
	return false;
};

bool
ExtendSignificance::checkMergedGroups(const std::string &t_prefix, Aux::Timer &timer,
                                      const GroupStubs &groupStubs, index statPosInt, node v,
                                      std::vector<std::pair<double, node>> &groupEdges,
                                      const std::vector<double> &groupSigs, count calcedGroups,
                                      std::vector<std::pair<node, double>> &nodeScores) const {
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

	for (++it; it < groupEdges.begin() + calcedGroups; ++it) {
		node group = it->second;
		count newInternalStubs = 0;
		for (node w : mergedGroups) {
			newInternalStubs += 2 * (int) coarseGraph.weight(w, group);
		}
		// Groups have to be connected (?)
		if (newInternalStubs == 0)
			continue;
		// Merge group
		extStubsMerged -= groupStubs.groupTotal[group];
		totalStubsMerged += groupStubs.groupTotal[group];
		outStubsMerged += groupStubs.groupOutgoing[group] - newInternalStubs;
		extNodesMerged -= coarseSizes[group];
		numEdgesMerged += (int) it->first;
		mergedGroups.push_back(group);
//			addTime(timer, t_prefix + "g    Merge groups");

		// Calculate new significance
		double significance = calcScore(nodeDegree, numEdgesMerged, outStubsMerged,
		                                extStubsMerged, extNodesMerged, statPosInt);
		added = addIfSignificant(nodeScores, v, significance);
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
ExtendSignificance::calcScore(int nodeDegree, int kIn, int grOut, int extStubs, int extNodes,
                              int position) const {
	// TODO: bootInterval anschauen
	assert(grOut <= extStubs);
	if (subtractNodeDegree)
		extStubs -= nodeDegree; // Remove node stubs from the free external stubs
	else
		extStubs -= 1;

	double rScore = Stochastics::compute_simple_fitness(kIn, grOut, extStubs, nodeDegree);

	double returnVal = 1.0;
	if (useSigMemo && position == 0) {
		double minR = sigTable.getValue(extNodes);
		if (rScore <= minR)
			returnVal = rScore / 16;
	} else {
		returnVal = Stochastics::order_statistics_left_cumulative(
				extNodes, extNodes - position, rScore);

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