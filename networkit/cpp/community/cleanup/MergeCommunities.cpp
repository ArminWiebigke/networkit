/*
 * MergeCommunities.cpp
 *
 * Created: 2019-09-26
 * Author: Armin Wiebigke
 */

#include "MergeCommunities.h"

namespace NetworKit {

using Community = MergeCommunities::Community;

MergeCommunities::MergeCommunities(const Graph &graph, std::set<Community> discardedCommunities,
                                   SingleCommunityCleanUp &singleCommunityCleanUp)
		: graph(graph),
		  discardedCommunities(std::move(discardedCommunities)),
		  stochastic(graph.numberOfNodes() + 2 * graph.numberOfEdges()),
		  singleCommunityCleanUp(singleCommunityCleanUp) {
}

void MergeCommunities::run() {
	int iterations = 2;
	for (count i = 0; i < iterations; ++i) {
		createDiscardedCommunitiesGraph();
		mergeCommunities();
		checkMergedCommunities();
	}
	hasRun = true;
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
	outgoingGroupStubs.clear();
	outgoingGroupStubs.resize(numDiscardedCommunities);
	totalGroupStubs.clear();
	totalGroupStubs.resize(numDiscardedCommunities);
	totalStubs = 0;
	graph.forEdges([&](node u, node v, edgeweight weight) {
		auto comms1 = nodeMemberships[u];
		auto comms2 = nodeMemberships[v];
		for (index comm1 : comms1) {
			for (index comm2 : comms2) {
				discardedCommunitiesGraph.increaseWeight(comm1, comm2, weight);
				totalStubs += 2;
				++totalGroupStubs[comm1];
				++totalGroupStubs[comm2];
				if (comm1 != comm2) {
					++outgoingGroupStubs[comm1];
					++outgoingGroupStubs[comm2];
				}
			}
		}
	});
}

void MergeCommunities::mergeCommunities() {
	mergedCommunities = Partition(discardedCommunitiesGraph.numberOfNodes());
	mergedCommunities.allToSingletons();
	count maxIterations = 20;
	for (count i = 0; i < maxIterations; ++i) {
		INFO("Merge communities iteration ", i);
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
	std::map<index, double> edgesToCommunities;
	discardedCommunitiesGraph.forEdgesOf(u, [&](node, node v, edgeweight weight) {
		auto neighborCommunity = mergedCommunities.subsetOf(v);
		edgesToCommunities[neighborCommunity] += weight;
	});

	index currentCommunity = mergedCommunities.subsetOf(u);
	edgeweight loop = discardedCommunitiesGraph.weight(u, u);
	count degree = discardedCommunitiesGraph.weightedDegree(u);
	count degreeWithoutLoop = degree - 2 * loop;

	index bestNeighborCommunity = currentCommunity;
	double bestScore = 1.;
	count stubsIntoCurrent = 0;
	count stubsIntoNew = 0;
	for (const auto &edgesToCommunity : edgesToCommunities) {
		count numEdgesIntoCommunity;
		index neighborCommunity;
		std::tie(neighborCommunity, numEdgesIntoCommunity) = edgesToCommunity;
		double score = 1.;
		count externalStubs = totalStubs - totalGroupStubs[neighborCommunity];
		if (neighborCommunity == currentCommunity) {
			// Calculate s-Score as if the node was not in the community
			numEdgesIntoCommunity -= loop;
			stubsIntoCurrent = numEdgesIntoCommunity;
			count outgoingFromNode = degree - 2 * loop - numEdgesIntoCommunity;
			count commOut =
					outgoingGroupStubs[currentCommunity] - outgoingFromNode + numEdgesIntoCommunity;
			externalStubs += degree;
			bool communityIsEmpty = commOut != 0;
			if (!communityIsEmpty) {
				score = stochastic.rScore(degree, numEdgesIntoCommunity,
				                          commOut, externalStubs);
			}
		} else {
			score = stochastic.rScore(degree, numEdgesIntoCommunity,
			                          outgoingGroupStubs[neighborCommunity],
			                          externalStubs);
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
		outgoingGroupStubs[currentCommunity] += -degreeWithoutLoop + 2 * stubsIntoCurrent;
		outgoingGroupStubs[bestNeighborCommunity] += degreeWithoutLoop - 2 * stubsIntoNew;
		totalGroupStubs[currentCommunity] -= degree;
		totalGroupStubs[bestNeighborCommunity] += degree;
		return true;
	}
	return false;
}

void MergeCommunities::checkMergedCommunities() {
	index communityId = 0;
	for (const auto &communitiesToMerge : mergedCommunities.getSubsets()) {
		if (communityId++ % (mergedCommunities.numberOfSubsets() / 10) == 0) {
			INFO("Clean merged community ", communityId, "/", mergedCommunities.numberOfSubsets());
		}
		if (communitiesToMerge.size() == 1)
			continue;
		Community mergedCommunity;
		for (index communityId : communitiesToMerge) {
			auto community = coarseToFineMapping[communityId];
			discardedCommunities.erase(community);
			for (node u : community)
				mergedCommunity.insert(u);
		}
		Community cleanedCommunity = singleCommunityCleanUp.clean(mergedCommunity);
		if (cleanedCommunity.empty())
			discardedCommunities.insert(mergedCommunity);
		else
			cleanedCommunities.push_back(cleanedCommunity);
	}
}

std::vector<Community> MergeCommunities::getCleanedCommunities() {
	return cleanedCommunities;
}

std::string MergeCommunities::toString() const {
	return "MergeCommunities";
}

bool MergeCommunities::isParallel() const {
	return false;
}

} /* namespace NetworKit */
