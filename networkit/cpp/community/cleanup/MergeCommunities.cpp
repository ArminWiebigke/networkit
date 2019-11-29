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
				   StochasticDistribution& stochasticDistribution,
                                   double significanceThreshold,
                                   double scoreThreshold,
                                   double minOverlapRatio,
                                   count maxCommunitySize)
		: graph(graph),
		  discardedCommunities(std::move(discardedCommunities)),
		  stochasticDistribution(stochasticDistribution),
		  stochastic(stochasticDistribution),
		  significanceThreshold(significanceThreshold),
		  scoreThreshold(scoreThreshold),
		  minOverlapRatio(minOverlapRatio),
		  maxCommunitySize(maxCommunitySize) {
}

void MergeCommunities::run() {
	int iterations = 2;
	for (count i = 0; i < iterations; ++i) {
		createDiscardedCommunitiesGraph();
		tryToMergeCommunities();
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

	graph.forNodes([&](node u) {
		const auto& comms1 = nodeMemberships[u];
		if (comms1.empty()) return;

		graph.forNeighborsOf(u, [&](node, node v, edgeweight weight) {
			if (u > v) return;
			const auto& comms2 = nodeMemberships[v];
			for (index comm2 : comms2) { // comms2 might be empty, but we have already checked comms1
				for (index comm1 : comms1) {
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
	});

	if (totalStubs + numDiscardedCommunities > stochasticDistribution.maxValue()) {
		stochasticDistribution.setMaxValue(totalStubs + numDiscardedCommunities);
	}
}

void MergeCommunities::tryToMergeCommunities() {
	mergedCommunities = Partition(discardedCommunitiesGraph.numberOfNodes());
	mergedCommunities.allToSingletons();
	count maxIterations = 20;
	for (count i = 0; i < maxIterations; ++i) {
		DEBUG("Merge communities iteration ", i);
		count nodesChanged = 0;
		discardedCommunitiesGraph.forNodesInRandomOrder([&](node u) {
			bool wasMoved = tryLocalMove(u);
			if (wasMoved)
				++nodesChanged;
		});
		if (nodesChanged == 0)
			break;
	}
	INFO("Merged ", mergedCommunities.numberOfElements(), " discarded communities into ",
	     mergedCommunities.numberOfSubsets(), " communities");
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
			// Calculate r-Score as if the node was not in the community
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
	// Store number of communities as this is not an O(1) lookup
	const count numMergedCommunities = mergedCommunities.numberOfSubsets();
	INFO("Check significance of ", numMergedCommunities, " merged communities");
	count skippedCommunities = 0;

	std::vector<index> partitionMap(mergedCommunities.upperBound(), none);
	std::vector<std::vector<node>> mergedCommunitiesSubsets;
	mergedCommunitiesSubsets.reserve(numMergedCommunities);
	mergedCommunities.forEntries([&](index e, index s){
		if (partitionMap[s] == none) {
			partitionMap[s] = mergedCommunitiesSubsets.size();
			mergedCommunitiesSubsets.emplace_back();
		}
		mergedCommunitiesSubsets[partitionMap[s]].push_back(e);
	});

	#pragma omp parallel
	{
		SingleCommunityCleanUp singleCommunityCleanUp(graph, stochasticDistribution, scoreThreshold, significanceThreshold, minOverlapRatio);
#pragma omp for schedule(dynamic, 1)
		for (omp_index i = 0; i < static_cast<omp_index>(mergedCommunitiesSubsets.size()); ++i) {
			const std::vector<node>& communitiesToMerge(mergedCommunitiesSubsets[i]);

			DEBUG("Clean merged community ", ++communityCount, "/", numMergedCommunities);
			if (communitiesToMerge.size() == 1)
				continue;

			Community mergedCommunity;

			for (index communityId : communitiesToMerge) {
				const auto& community = coarseToFineMapping[communityId];
#pragma omp critical
				discardedCommunities.erase(community);
				for (node u : community) {
					mergedCommunity.insert(u);
				}
				if (mergedCommunity.size() >= maxCommunitySize)
					break;
			}
			if (mergedCommunity.size() >= maxCommunitySize) {
#pragma omp atomic update
				++skippedCommunities;
				continue;
			}
			Community cleanedCommunity = singleCommunityCleanUp.clean(mergedCommunity);
#pragma omp critical
			{
				if (cleanedCommunity.empty())
					discardedCommunities.emplace(std::move(mergedCommunity));
				else
					cleanedCommunities.emplace_back(std::move(cleanedCommunity));
			}
		}
	}
	INFO("Skipped ", skippedCommunities, " large communities (max size ", maxCommunitySize, ")");
}

std::vector<Community> MergeCommunities::getCleanedCommunities() {
	return cleanedCommunities;
}

std::string MergeCommunities::toString() const {
	return "MergeCommunities";
}

bool MergeCommunities::isParallel() const {
	return true;
}

} /* namespace NetworKit */
