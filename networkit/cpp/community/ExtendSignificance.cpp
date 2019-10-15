/*
 * ExtendSignificance.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#include <set>

#include "ExtendSignificance.h"
#include "../auxiliary/ParseString.h"
#include "../coarsening/ParallelPartitionCoarsening.h"
#include "../oslom/Stochastics.h"
#include "../cleanup/StochasticDistribution.h"
#include "EgoSplitting.h"

#define W(x) #x << "=" << x << ", "
#define true_or_throw(cond, msg) if (!cond) throw std::runtime_error(msg)

namespace NetworKit {

ExtendSignificance::ExtendSignificance(const EgoNetData &egoNetData,
                                       const Partition &basePartition, count maxCandidates,
                                       const Graph &egoGraph, node egoNode)
		: ExtendEgoNetStrategy(egoNetData, maxCandidates, egoGraph, egoNode),
		  basePartition(basePartition),
		  sigTable(egoNetData.sigTable),
		  stochasticSignificance(egoNetData.stochasticSignificance),
		  useSigMemo(parameters.at("useSigMemo") == "Yes"),
		  mergeGroups(parameters.at("signMerge") == "Yes"),
		  sortGroupsStrat(parameters.at("sortGroups") == "Significance"),
		  maxSignificance(Aux::stringToDouble(parameters.at("maxSignificance"))),
		  maxGroupCnt(std::stoi(parameters.at("maxGroupsConsider"))),
		  minEdgesToGroup(std::stoi(parameters.at("minEdgesToGroupSig"))),
		  significantGroup(egoNetData.significantGroup),
		  edgesToGroups(egoNetData.edgesToGroups),
		  onlyCheckMaxCandidates(parameters.at("onlyCheckSignOfMaxCandidates") == "Yes") {
	if (basePartition.numberOfSubsets() == 0)
		throw std::runtime_error("Missing Base Partition for Significance!");
	setMemoizationFunction();
}

void ExtendSignificance::setMemoizationFunction() const {
	double maxSig = maxSignificance;
	// For a given number of external nodes, calculate the minimal s-score that is significant
	// We can then simply compare the s-values instead of calculating the ordered statistics
	// Only works if the position is known before. 0 ( = best random node) in this case.
	const StochasticSignificance &stoch = stochasticSignificance;
	sigTable.trySetValueFunc([maxSig, &stoch](index extNodes) {
		double minRScore = 0.125;
		double step = 0.0625;
		count iterations = 20;
		for (count i = 0; i < iterations; ++i) {
//			double sig = Stochastics::order_statistics_left_cumulative(extNodes, extNodes,
//			                                                           minRScore);
			double sig = stoch.orderStatistic(minRScore, extNodes, 1);
			if (sig < maxSig) {
				minRScore += step;
			} else {
				minRScore -= step;
			}
			step /= 2;
		}
		return minRScore;
	});
}

void ExtendSignificance::run() {
	Aux::Timer timer;
	timer.start();
//	Stochastics::init(G.numberOfEdges());

	// Get coarse graph from base partition
	createCoarseGraph();
	addTime(timer, "0    Coarsening");

	std::vector<node> candidates = getCandidates();
	addTime(timer, "3    Find candidates");

	addExtEdges(candidates); // TODO: Is this only for checking/debugging?
	groupProperties = calcGroupStubsCounts();
	addTime(timer, "4    Calc outgoing edges");

	candidatesSorted = sortCandidatesByEdges(candidates);
	addTime(timer, "5    Sort candidates");

	if (onlyCheckMaxCandidates) {
		double evalFactor = Aux::stringToDouble(parameters.at("Check Candidates Factor"));
		count evalCandidates = evalFactor * maxExtendedNodes;
		if (candidatesSorted.size() > evalCandidates)
			candidatesSorted.resize(evalCandidates);
	}

	auto debugTest = [&]() {
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
	assert(debugTest());

	findSignificantCandidates("6");
	addTime(timer, "6    Calc significance one");

	count iterations = std::stoi(parameters.at("secondarySigExtRounds"));
	for (count i = 0; i < iterations; ++i) {
		secondRound();
	}
	addTime(timer, "8    Check updated candidates");

#ifndef NDEBUG
	std::set<node> resultSet;
	for (node v : significantCandidates) {
		assert(resultSet.count(v) == 0);
		resultSet.insert(v);
	}
	assert(resultSet.size() == significantCandidates.size());
#endif

	resetData();
	addTime(timer, "9    Reset Data");

	hasRun = true;
}

void ExtendSignificance::updateCandidates() {
	std::unordered_set<node> updatedCandidates;
	for (node w : addedCandidates) {
		node group = significantGroup[w];
		// Update group stubs
		groupProperties.groupTotal[group] += G.degree(w);
		groupProperties.groupOutgoing[group] += G.degree(w) - 2 * edgesToGroups[w][group];
		groupProperties.externalNodes[group] -= 1;
		groupProperties.externalStubs[group] -= G.degree(w);
		G.forNeighborsOf(w, [&](node v) {
			// Update other candidates
			if (!edgesToGroups[v].empty()) {
				// v is a neighbor of neighbor
				edgesToGroups[v][group] += 1;
				updatedCandidates.insert(v);
			}
			if (significantGroup[v] == group) {
				// v is in the same group
				count newInternalStubs = 2;
				if (addedCandidates.count(v)) {
					// Don't count this edge twice if we add both candidates at the same time
					newInternalStubs = 1;
				}
				groupProperties.groupOutgoing[group] -= newInternalStubs;
			}
		});
		edgesToGroups[w].clear();
	}
	bool onlyUseUpdatedCandidates = parameters.at("onlyUpdatedCandidates") == "Yes";
	auto newEnd = std::remove_if(
			candidatesSorted.begin(), candidatesSorted.end(),
			[&](node candidate) {
				return addedCandidates.count(candidate) > 0
				       || (onlyUseUpdatedCandidates && updatedCandidates.count(candidate) == 0);
			});
	candidatesSorted.resize(newEnd - candidatesSorted.begin());
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
	if (addedCandidates.empty())
		return;
	updateCandidates();
	addedCandidates.clear();
	findSignificantCandidates("8"); // TODO: Remove timer?
}

/**
 * Sort candidates by number of edges into the egoNet, filter bad candidates (less than 3 edges)
 */
std::vector<node>
ExtendSignificance::sortCandidatesByEdges(const std::vector<node> &candidates) const {
	std::vector<std::pair<count, node>> candidatesAndEdgeCnt;
	candidatesAndEdgeCnt.reserve(candidates.size());
	count minEdges = 3;
	for (node w : candidates) {
		count edgesIntoEgo = std::accumulate(edgesToGroups[w].begin(), edgesToGroups[w].end(),
		                                     (count) 0);
		if (edgesIntoEgo >= minEdges)
			candidatesAndEdgeCnt.emplace_back(edgesIntoEgo, w);
	}
	std::sort(candidatesAndEdgeCnt.rbegin(), candidatesAndEdgeCnt.rend());

	std::vector<node> sortedCandidates(candidatesAndEdgeCnt.size());
	for (index i = 0; i < candidatesAndEdgeCnt.size(); ++i) {
		sortedCandidates[i] = candidatesAndEdgeCnt[i].second;
	}
	return sortedCandidates;
}

/**
 * Insert a node that represents all external nodes and add the outgoing edges
 */
node
ExtendSignificance::addExtEdges(std::vector<node> &candidates) {
	node externalNode = coarseGraph.addNode();
	assert(numGroups == externalNode);
	coarseSizes.push_back(G.numberOfNodes() - egoGraph.numberOfNodes());
	// Calc outgoing edges for each group
	candidates.push_back(egoNode); // Add edges to egoNode, but remove it as candidate later
	for (node w : candidates) {
		for (index group = 0; group < edgesToGroups[w].size(); ++group) {
			count numEdges = edgesToGroups[w][group];
			coarseGraph.increaseWeight(group, externalNode, numEdges);
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
			edgesToGroups.insert(w, std::vector<count>(numGroups));
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

	for (node v : directNeighbors) {
		node vLoc = egoToCoarse[egoMapping.local(v)];
		G.forEdgesOf(v, [&](node, node w, edgeweight weight) {
			countEdge(vLoc, w, weight);
		});
	}
	removeNeighborCandidates();
	removeEgoNodeCandidate(candidates);
	return candidates;
}

void ExtendSignificance::removeEgoNodeCandidate(std::vector<node> &candidates) {
	for (index i = 0; i < candidates.size(); ++i) {
		if (candidates[i] == egoNode) {
			candidates[i] = candidates.back();
			candidates.pop_back();
			break;
		}
	}
	edgesToGroups[egoNode].clear();
}

void
ExtendSignificance::findSignificantCandidates(const std::string &t_prefix) {
	Aux::Timer timer;
	timer.start();

	for (node candidate : candidatesSorted) {
		if (enoughSignificantCandidates())
			break;
		checkCandidate(t_prefix, timer, candidate);
	}
}

bool ExtendSignificance::enoughSignificantCandidates() const {
	return significantCandidates.size() >= maxExtendedNodes;
}

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

/**
 * For each group, calculate its properties (internal/outgoing stubs, external stubs/nodes)
 * @return data structure that contains the properties of all groups
 */
GroupProperties ExtendSignificance::calcGroupStubsCounts() const {
	GroupProperties stubs(numGroups);
	for (index group = 0; group < numGroups; ++group) {
		if (!coarseGraph.hasNode(group))
			continue;
		auto groupDegree = (count) coarseGraph.weightedDegree(group);
		stubs.groupTotal[group] = groupDegree;
		stubs.groupOutgoing[group] = groupDegree - 2 * (count) coarseGraph.weight(group, group);
		stubs.externalStubs[group] = 2 * G.numberOfEdges() - groupDegree;
		stubs.externalNodes[group] = G.numberOfNodes() - coarseSizes[group];
	}
	return stubs;
}

/**
 * Check the significance of the candidate to all groups. If the candidate is significant,
 * add it to the list of significant candidates.
 * @param v
 * @param groupEdges
 * @param groupSigs
 * @return true if the candidate is significant to at least on of the groups, else false
 */
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
		count numEdges = it.first;
		node p = it.second;
		double significance = calcScore(nodeDegree, numEdges, groupProperties.groupOutgoing[p],
		                                groupProperties.externalStubs[p],
		                                groupProperties.externalNodes[p]);
		groupSigs[p] = significance;
		added = addIfSignificant(v, significance, p);
		if (added)
			break;
	}
	return added;
}

bool ExtendSignificance::addIfSignificant(node v, double significance, node group) {
	if (significance <= maxSignificance) {
		significantCandidates.push_back(v);
		significantGroup.insert(v, group);
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
	addTime(timer, t_prefix + "d    Sort groups by Significance");

	auto it = groupEdges.begin();
	node bestGroup = it->second;
	std::vector<node> mergedGroups{bestGroup};
	count totalStubsMerged = groupProperties.groupTotal[bestGroup];
	count outStubsMerged = groupProperties.groupOutgoing[bestGroup];
	count extStubsMerged = groupProperties.externalStubs[bestGroup];
	count extNodesMerged = groupProperties.externalNodes[bestGroup];
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
		extStubsMerged -= groupProperties.groupTotal[group];
		totalStubsMerged += groupProperties.groupTotal[group];
		outStubsMerged += groupProperties.groupOutgoing[group] - newInternalStubs;
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
ExtendSignificance::calcScore(count nodeDegree, count kIn, count grOut, count groupExtStubs,
                              count extNodes) const {
//	double rScore = Stochastics::compute_simple_fitness(kIn, grOut, openStubs, nodeDegree);
	double rScore = stochasticSignificance.rScore(nodeDegree, kIn, grOut, groupExtStubs);
	double returnVal = 1.0;
	if (useSigMemo) {
		// TODO: Why does this add more nodes?
		double minRScore = sigTable.getValue(extNodes);
		if (rScore <= minRScore) {
			returnVal = rScore;
			assert(stochasticSignificance.orderStatistic(rScore, extNodes, 1) <=
			       1.01 * maxSignificance);
//			assert(Stochastics::order_statistics_left_cumulative(
//					extNodes, extNodes - orderStatPos, rScore) <= maxSignificance);
		} else {
			assert(stochasticSignificance.orderStatistic(rScore, extNodes, 1) >=
			       0.99 * maxSignificance);
//			assert(Stochastics::order_statistics_left_cumulative(
//					extNodes, extNodes - orderStatPos, rScore) >= maxSignificance);
		}
	} else {
		returnVal = stochasticSignificance.orderStatistic(rScore, extNodes, 1);
//		returnVal = Stochastics::order_statistics_left_cumulative(
//				extNodes, extNodes - orderStatPos, rScore);
	}
	return returnVal;
}

std::string ExtendSignificance::toString() const {
	return "ExtendSignificance";
}

bool ExtendSignificance::isParallel() const {
	return false;
}

void ExtendSignificance::resetData() {
	significantGroup.reset();
	edgesToGroups.reset();
}
} /* namespace NetworKit */