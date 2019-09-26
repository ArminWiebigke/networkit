/*
 * MergeCommunities.cpp
 *
 * Created: 2019-09-26
 * Author: Armin Wiebigke
 */

#include "MergeCommunities.h"

namespace NetworKit {

using Community = MergeCommunities::Community;

std::vector<Community> & MergeCommunities::getCleanedCommunities() {
	return cleanedCommunities;
}

void MergeCommunities::run() {
	for (count i = 0; i < 2; ++i) {
		createDiscardedCommunitiesGraph();
		localOptimization();
		checkMergedCommunities();
	}
}

void MergeCommunities::createDiscardedCommunitiesGraph() {
	count numDiscardedCommunities = discardedCommunities.size();
	discardedCommunitiesGraph = Graph(numDiscardedCommunities, true, false);

	// Get node memberships
	std::vector<std::vector<index>> nodeMemberships(graph.upperNodeIdBound());
	index commId = 0;
	coarseToFineMapping.resize(numDiscardedCommunities);
	for (const auto &community : discardedCommunities) {
		for (node u : community) {
			nodeMemberships[u].push_back(commId);
			coarseToFineMapping[commId].insert(u);
		}
		++commId;
	}

	// Count edges between communities
	cOut.clear();
	cOut.resize(numDiscardedCommunities);
	cTotal.clear();
	cTotal.resize(numDiscardedCommunities);
	totalStubs = 0;
	graph.forEdges([&](node u, node v, edgeweight weight) {
		auto comms1 = nodeMemberships[u];
		auto comms2 = nodeMemberships[v];
		for (index comm1 : comms1) {
			for (index comm2 : comms2) {
				discardedCommunitiesGraph.increaseWeight(comm1, comm2, weight);
				totalStubs += 2;
				++cTotal[comm1];
				++cTotal[comm2];
				if (comm1 != comm2) {
					++cOut[comm1];
					++cOut[comm2];
				}
			}
		}
	});
}

void MergeCommunities::localOptimization() {
	mergedCommunities = Partition(discardedCommunitiesGraph.numberOfNodes());
	mergedCommunities.allToSingletons();
	count maxIterations = 20;
	for (count i = 0; i < maxIterations; ++i) {
		count nodesChanged = 0;
		discardedCommunitiesGraph.forNodesInRandomOrder([&](node u) {
			bool wasMoved = tryLocalMove(u);
			if (wasMoved)
				++nodesChanged;
		});
		if (nodesChanged == 0)
			break;
	}
}

bool MergeCommunities::tryLocalMove(node u) {
	std::map<index, double> kIn;
	discardedCommunitiesGraph.forEdgesOf(u, [&](node, node v, edgeweight weight) {
		auto neighborCommunity = mergedCommunities.subsetOf(v);
		kIn[neighborCommunity] += weight;
	});

	index currentCommunity = mergedCommunities.subsetOf(u);
	edgeweight loop = discardedCommunitiesGraph.weight(u, u);
	count degree = discardedCommunitiesGraph.weightedDegree(u);
	count degreeWithoutLoop = degree - 2 * loop;

	index bestNeighborCommunity = currentCommunity;
	double bestScore = 2.;
	count stubsIntoCurrent = 0;
	count stubsIntoNew = 0;
	for (const auto &x : kIn) {
		count numEdgesIntoCommunity;
		index neighborCommunity;
		std::tie(neighborCommunity, numEdgesIntoCommunity) = x;
		double score;
		count extStubs = totalStubs - cTotal[neighborCommunity];
		if (neighborCommunity == currentCommunity) {
			numEdgesIntoCommunity -= loop;
			stubsIntoCurrent = numEdgesIntoCommunity;
			count outgoingFromThis = degree - 2 * loop - numEdgesIntoCommunity;
			count commOut = cOut[currentCommunity] - outgoingFromThis + numEdgesIntoCommunity;
			extStubs += degree;
			if (commOut == 0) {
				score = 1.;
			} else {
				auto scoreTuple = stochastic.sScore(degree, numEdgesIntoCommunity,
				                                    commOut, extStubs);
				score = scoreTuple.first;
			}
		} else {
			auto scoreTuple = stochastic.sScore(degree, numEdgesIntoCommunity,
			                                    cOut[neighborCommunity], extStubs);
			score = scoreTuple.first;
		}
		if (score < bestScore) {
			stubsIntoNew = numEdgesIntoCommunity;
			bestScore = score;
			bestNeighborCommunity = neighborCommunity;
		}
	}

	if (bestNeighborCommunity != currentCommunity) {
		// Move node into new community
		mergedCommunities.moveToSubset(bestNeighborCommunity, u);
		cOut[currentCommunity] += -degreeWithoutLoop + 2 * stubsIntoCurrent;
		cOut[bestNeighborCommunity] += degreeWithoutLoop - 2 * stubsIntoNew;
		cTotal[currentCommunity] -= degree;
		cTotal[bestNeighborCommunity] += degree;
		return true;
	}
	return false;
}

void MergeCommunities::checkMergedCommunities() {
	for (const auto &communities : mergedCommunities.getSubsets()) {
		if (communities.size() == 1)
			continue;
		Community mergedCommunity;
		for (auto community : communities) {
			auto nodes = coarseToFineMapping[community];
			discardedCommunities.erase(nodes);
			for (node u : nodes)
				mergedCommunity.insert(u);
		}
		Community cleanedCommunity = singleCommunityCleanUp.clean(mergedCommunity);
//		std::cout << "Merged: " << mergedCommunity.size() << " -> " << cleanedCommunity.size()
//		          << std::endl;
		if (cleanedCommunity.empty())
			discardedCommunities.insert(mergedCommunity);
		else
			cleanedCommunities.push_back(cleanedCommunity);
	}
}
} /* namespace NetworKit */
