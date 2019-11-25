/*
 * ExtendSignificance.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#include <set>

#include "ExtendSignificance.h"
#include "../../auxiliary/ParseString.h"
#include "../../coarsening/ParallelPartitionCoarsening.h"
#include "../../oslom/Stochastics.h"
#include "../cleanup/StochasticDistribution.h"
#include "EgoSplitting.h"

#define W(x) #x << "=" << x << ", "
#define true_or_throw(cond, msg) if (!cond) throw std::runtime_error(msg)

namespace NetworKit {

ExtendSignificance::ExtendSignificance(EgoNetData &egoNetData,
                                       const Partition &basePartition, count maxCandidates,
                                       const Graph &egoGraph, node egoNode)
		: ExtendEgoNetStrategy(egoNetData, maxCandidates, egoGraph, egoNode),
		  basePartition(basePartition),
		  sigTable(egoNetData.sigTable),
		  stochasticSignificance(egoNetData.stochasticSignificance),
		  useSigMemo(parameters.at("useSigMemo") == "Yes"),
		  mergeGroups(parameters.at("signMerge") == "Yes"),
		  sortGroupsBySignificance(parameters.at("sortGroups") == "Significance"),
		  maxSignificance(Aux::stringToDouble(parameters.at("maxSignificance"))),
		  maxGroupCnt(std::stoi(parameters.at("maxGroupsConsider"))),
		  minEdgesToGroup(std::stoi(parameters.at("minEdgesToGroupSig"))),
		  significantGroup(egoNetData.significantGroup),
		  edgesToGroups(egoNetData.edgesToGroups),
		  onlyCheckMaxCandidates(parameters.at("onlyCheckSignOfMaxCandidates") == "Yes") {
	if (basePartition.numberOfSubsets() == 0 && egoGraph.numberOfNodes() > 0)
		throw std::runtime_error("Missing Base Partition for Significance!");
	if (significantGroup.upperBound() < G.upperNodeIdBound())
		significantGroup.setUpperBound(G.upperNodeIdBound());
	if (edgesToGroups.upperBound() < G.upperNodeIdBound())
		edgesToGroups.setUpperBound(G.upperNodeIdBound());
	setMemoizationFunction();
}

void ExtendSignificance::setMemoizationFunction() const {
	if (!sigTable.valueFunctionIsSet()) {
		double maxSig = maxSignificance;
		// For a given number of external nodes, calculate the minimal s-score that is significant
		// We can then simply compare the s-values instead of calculating the ordered statistics
		// Only works if the position is known before. 1 ( = best random node) in this case.
		const StochasticSignificance &stoch = stochasticSignificance;
		sigTable.setValueFunc([maxSig, &stoch](index extNodes) {
			double minRScore = 0.125;
			double step = 0.0625;
			count iterations = 20;
			for (count i = 0; i < iterations; ++i) {
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
}

void ExtendSignificance::run() {
	Aux::Timer timer;
	timer.start();

	createCoarseGraph();
	addTime(timer, "0    Coarsening");

	auto candidates = getCandidates();
	addTime(timer, "3    Find candidates");

	insertOutgoingEdgesIntoCoarseGraph(candidates);
	groupProperties = calculateGroupProperties();
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

	count updateIterations = std::stoi(parameters.at("secondarySigExtRounds"));
	for (count i = 0; i < updateIterations; ++i) {
		if (addedCandidates.empty())
			break;
		updateCandidates();
		addedCandidates.clear();
		findSignificantCandidates("6");
	}
	addTime(timer, "6    Find significant candidates");

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

/**
 * Remove significant nodes as candidates and update the remaining candidates
 */
void ExtendSignificance::updateCandidates() {
	std::unordered_set<node> updatedCandidates;
	for (node w : addedCandidates) {
		node group = significantGroup[w];
		// Update group stubs
		GroupProperties &properties = groupProperties[group];
		properties.groupTotal += G.degree(w);
		int oldOutgoing = properties.groupOutgoing;
		assert(oldOutgoing + (int) G.degree(w) - 2 * (int) edgesToGroups[w][group] >= 0);
		properties.groupOutgoing += G.degree(w) - 2 * edgesToGroups[w][group];
		assert(properties.groupOutgoing < 2 * G.numberOfEdges());
		properties.externalNodes -= 1;
		properties.externalStubs -= G.degree(w);
		G.forNeighborsOf(w, [&](node v) {
			// Update other candidates
			if (!edgesToGroups[v].empty()) {
				// v is a candidate
				edgesToGroups[v][group] += 1;
				updatedCandidates.insert(v);
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

/**
 * Get coarse graph from base partition
 */
void ExtendSignificance::createCoarseGraph() {
	ParallelPartitionCoarsening coarsening(egoGraph, basePartition, true, false);
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
 * Insert a node in the coarse graph that represents all external nodes and add the outgoing edges
 * of each group.
 */
void
ExtendSignificance::insertOutgoingEdgesIntoCoarseGraph(std::vector<node> &candidates) {
	std::vector<count> groupToExternal(coarseGraph.numberOfNodes(), 0);
	auto countEdgesOfNode = [&](node v) {
		for (index group = 0; group < edgesToGroups[v].size(); ++group) {
			count numEdges = edgesToGroups[v][group];
			groupToExternal[group] += numEdges;
		}
	};
	countEdgesOfNode(egoNode);
	for (node v : candidates) {
		countEdgesOfNode(v);
	}
	node externalNode = coarseGraph.addNode();
	coarseSizes.push_back(G.numberOfNodes() - egoGraph.numberOfNodes());
	for (index i = 0; i < coarseGraph.numberOfNodes() - 1; ++i) {
		coarseGraph.addEdge(i, externalNode, groupToExternal[i]);
	}
}

std::vector<node>
ExtendSignificance::getCandidates() {
	std::vector<node> candidates;
	auto countEdge = [&](node group, node candidate, edgeweight weight) {
		if (!edgesToGroups.indexIsUsed(candidate)) {
			candidates.push_back(candidate);
			edgesToGroups.insert(candidate, std::vector<count>(numGroups, 0));
		}
		edgesToGroups[candidate][group] += (count) weight;
	};
	std::vector<node> egoNetNodes = egoMapping.globalNodes();
	for (node egoNetNode : egoNetNodes) {
		node group = egoToCoarse[egoMapping.toLocal(egoNetNode)];
		G.forEdgesOf(egoNetNode, [&](node, node candidate, edgeweight weight) {
			countEdge(group, candidate, weight);
		});
	}

	auto shouldBeRemoved = [&](node x) {
		// Remove ego-node and nodes of the ego-net as extension candidates
		return egoMapping.isMapped(x) || x == egoNode;
	};
	auto removeEgoNetNodesAsCandidates = [&]() {
		for (node v : candidates) {
			if (shouldBeRemoved(v)) {
				edgesToGroups[v].clear();
			}
		}
		auto newEnd = std::remove_if(candidates.begin(), candidates.end(), [&](node w) {
			return (shouldBeRemoved(w));
		});
		candidates.resize(newEnd - candidates.begin());
	};
	removeEgoNetNodesAsCandidates();
	return candidates;
}

void
ExtendSignificance::findSignificantCandidates(const std::string &t_prefix) {
	for (node candidate : candidatesSorted) {
		if (enoughSignificantCandidates())
			break;
		checkCandidate(t_prefix, candidate);
	}
}

bool ExtendSignificance::enoughSignificantCandidates() const {
	return significantCandidates.size() >= maxExtendedNodes;
}

void
ExtendSignificance::checkCandidate(const std::string &t_prefix, node candidate) {
	Aux::Timer timer{};
	timer.start();

	std::vector<std::pair<double, node>> numEdgesToGroups = sortGroupsByEdges(candidate);
	if (numEdgesToGroups.empty())
		return;
	addTime(timer, t_prefix + "b    Sort groups by edge count");

	std::vector<double> significanceToGroups(numGroups);
	bool added = checkSignificanceToSingleGroups(candidate, numEdgesToGroups, significanceToGroups);
	addTime(timer, t_prefix + "c    Significance single");
	if (added)
		return;

	if (mergeGroups) {
		checkSignificanceToMergedGroups(t_prefix, timer, candidate, numEdgesToGroups,
		                                significanceToGroups);
		addTime(timer, t_prefix + "f    Significance merge");
	}
}

std::vector<std::pair<double, node>> ExtendSignificance::sortGroupsByEdges(node candidate) const {
	std::vector<std::pair<double, node>> groupEdges;
	groupEdges.reserve(edgesToGroups[candidate].size());
	for (node group = 0; group < edgesToGroups[candidate].size(); ++group) {
		count numEdges = edgesToGroups[candidate][group];
		if (numEdges >= minEdgesToGroup)
			groupEdges.emplace_back(numEdges, group);
	}
	std::sort(groupEdges.rbegin(), groupEdges.rend());
	return groupEdges;
}

/**
 * For each group, calculate its properties (internal/outgoing stubs, external stubs/nodes)
 * @return data structure that contains the properties of all groups
 */
std::vector<GroupProperties> ExtendSignificance::calculateGroupProperties() const {
	std::vector<GroupProperties> newGroupProperties(numGroups);
	for (index group = 0; group < numGroups; ++group) {
		if (!coarseGraph.hasNode(group))
			continue;
		auto groupDegree = (count) coarseGraph.weightedDegree(group);
		GroupProperties &properties = newGroupProperties[group];
		properties.groupTotal = groupDegree;
		properties.groupOutgoing = groupDegree - 2 * (count) coarseGraph.weight(group, group);
		assert(properties.groupOutgoing < 2 * G.numberOfEdges());
		properties.externalStubs = 2 * G.numberOfEdges() - groupDegree;
		properties.externalNodes = G.numberOfNodes() - coarseSizes[group];
	}
	return newGroupProperties;
}

/**
 * Check the significance of the candidate to all groups. If the candidate is significant,
 * add it to the list of significant candidates.
 * @param candidate
 * @param numEdgesToGroups
 * @param significanceToGroups significance of the candidate to the groups
 * @return true if the candidate is significant to at least on of the groups, else false
 */
bool ExtendSignificance::checkSignificanceToSingleGroups(
		node candidate,
		const std::vector<std::pair<double, node>> &numEdgesToGroups,
		std::vector<double> &significanceToGroups) {
	count nodeDegree = G.degree(candidate);
	bool added = false;
	count calcedGroups = 0;
	for (auto it : numEdgesToGroups) {
		if (calcedGroups >= maxGroupCnt)
			break;
		++calcedGroups;
		count numEdgesToGroup = it.first;
		node group = it.second;
		GroupProperties &properties = groupProperties[group];
		double sScore = calculateSScore(nodeDegree,
		                                numEdgesToGroup,
		                                properties.groupOutgoing,
		                                properties.externalStubs,
		                                properties.externalNodes);
		significanceToGroups[group] = sScore;
		added = addIfSignificant(candidate, sScore, group);
		if (added)
			break;
	}
	return added;
}

bool
ExtendSignificance::checkSignificanceToMergedGroups(
		const std::string &t_prefix, Aux::Timer &timer,
		node candidate,
		std::vector<std::pair<double, node>> &numEdgesToGroups,
		const std::vector<double> &groupSigs) {
	count nodeDegree = G.degree(candidate);
	if (sortGroupsBySignificance) {
		std::sort(numEdgesToGroups.begin(), numEdgesToGroups.end(),
		          [&](std::pair<double, node> a, std::pair<double, node> b) {
			          return (groupSigs[a.second]) > groupSigs[b.second];
		          });
		addTime(timer, t_prefix + "d    Sort groups by Significance");
	}

	auto it = numEdgesToGroups.begin();
	node bestGroup = it->second;
	std::vector<node> mergedGroups{bestGroup};
	GroupProperties &properties = groupProperties[bestGroup];
	count totalStubsMerged = properties.groupTotal;
	count outStubsMerged = properties.groupOutgoing;
	count extStubsMerged = properties.externalStubs;
	count extNodesMerged = properties.externalNodes;
	count numEdgesMerged = it->first;
	bool added = false;

	for (++it; it < numEdgesToGroups.end(); ++it) {
		node groupToMerge = it->second;
		if (groupSigs[groupToMerge] == 1.0)
			break;
		// Merge group
		count newInternalStubs = 0;
		for (node mergedGroup : mergedGroups) {
			newInternalStubs += 2 * (count) coarseGraph.weight(mergedGroup, groupToMerge);
		}
		GroupProperties &mergeProperties = groupProperties[groupToMerge];
		extStubsMerged -= mergeProperties.groupTotal;
		totalStubsMerged += mergeProperties.groupTotal;
		outStubsMerged += mergeProperties.groupOutgoing - newInternalStubs;
		extNodesMerged -= coarseSizes[groupToMerge];
		numEdgesMerged += (count) it->first;
		mergedGroups.push_back(groupToMerge);

		// Calculate new significance
		double sScore = calculateSScore(nodeDegree, numEdgesMerged, outStubsMerged,
		                                extStubsMerged, extNodesMerged);
		added = addIfSignificant(candidate, sScore, bestGroup);
		if (added) {
			break;
		}
	}
	return added;
}

double
ExtendSignificance::calculateSScore(count nodeDegree, count kIn, count grOut, count groupExtStubs,
                                    count extNodes) const {
	double rScore = stochasticSignificance.rScore(nodeDegree, kIn, grOut, groupExtStubs);
	double returnVal;
	if (useSigMemo) {
		double minRScore = sigTable.getValue(extNodes);
#ifndef NDEBUG
		if (rScore <= minRScore) {
			assert(stochasticSignificance.orderStatistic(rScore, extNodes, 1) <=
			       1.01 * maxSignificance);
		} else {
			assert(stochasticSignificance.orderStatistic(rScore, extNodes, 1) >=
			       0.99 * maxSignificance);
		}
#endif
		if (rScore <= minRScore) {
			returnVal = rScore;
		} else {
			returnVal = 1.0;
		}
	} else {
		returnVal = stochasticSignificance.orderStatistic(rScore, extNodes, 1);
	}
	return returnVal;
}

bool ExtendSignificance::addIfSignificant(node v, double significance, node group) {
	if (significance <= maxSignificance) {
		significantCandidates.push_back(v);
		significantGroup.insert(v, group);
		addedCandidates.insert(v);
		return true;
	}
	return false;
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