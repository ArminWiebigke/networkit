/*
 * SignificanceCommunityCleanUp.cpp
 *
 * Created: 2019-09-11
 * Author: Armin Wiebigke
 */

#include "SignificanceCommunityCleanUp.h"
#include "MergeCommunities.h"

namespace NetworKit {

using Community = SignificanceCommunityCleanUp::Community;

SignificanceCommunityCleanUp::SignificanceCommunityCleanUp(const Graph &graph,
                                                           const Cover &cover,
							   const StochasticDistribution& distribution,
                                                           double significanceThreshold,
                                                           double scoreThreshold,
                                                           double minOverlapRatio,
                                                           bool mergeDiscarded)
		: graph(graph),
		  cover(cover),
		  significanceThreshold(significanceThreshold),
		  scoreThreshold(scoreThreshold),
		  minOverlapRatio(minOverlapRatio),
		  mergeDiscarded(mergeDiscarded),
		  maxCommunitySize(0),
		  stochasticDistribution(distribution) {
}

void SignificanceCommunityCleanUp::run() {
	hasRun = false;
	cleanedCommunities = Cover(graph.upperNodeIdBound());
	cleanAllCommunities();
	if (mergeDiscarded) {
		mergeDiscardedCommunities();
	}
	hasRun = true;
}

Cover SignificanceCommunityCleanUp::getCover() {
	if (!hasFinished())
		throw std::runtime_error("Run the algorithm first!");
	return cleanedCommunities;
}

std::string SignificanceCommunityCleanUp::toString() const {
	return "SignificanceCommunityCleanUp";
}

bool SignificanceCommunityCleanUp::isParallel() const {
	return true;
}

void SignificanceCommunityCleanUp::cleanAllCommunities() {
	auto inputCommunities = cover.getSubsets();
	INFO("Clean ", inputCommunities.size(), " communities");
	#pragma omp parallel
	{
		SingleCommunityCleanUp singleCommunityCleanup(graph, stochasticDistribution, scoreThreshold, significanceThreshold, minOverlapRatio);
#pragma omp for schedule(dynamic, 1)
		for (omp_index i = 0; i < static_cast<omp_index>(inputCommunities.size()); ++i) {
			const Community& inputCommunity = inputCommunities[i];
			DEBUG("Clean community ", i, "/", inputCommunities.size(), " with size ", inputCommunity.size());
			auto cleanedCommunity = singleCommunityCleanup.clean(inputCommunity);
			#pragma omp critical
			{
				if (cleanedCommunity.empty() && mergeDiscarded)
					discardedCommunities.insert(inputCommunity);
				else {
					cleanedCommunities.addSubset(cleanedCommunity);
					maxCommunitySize = std::max(maxCommunitySize, cleanedCommunity.size());
				}
			}
		}
	}
}


void SignificanceCommunityCleanUp::mergeDiscardedCommunities() {
	INFO("Try to merge ", discardedCommunities.size(), " discarded communities");
	MergeCommunities mergeCommunities(graph, std::move(discardedCommunities), stochasticDistribution, significanceThreshold, scoreThreshold, minOverlapRatio,
	                                  2 * maxCommunitySize);
	mergeCommunities.run();
	for (const auto &community : mergeCommunities.getCleanedCommunities()) {
		cleanedCommunities.addSubset(community);
	}
}

} /* namespace NetworKit */
