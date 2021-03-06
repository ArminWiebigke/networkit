/*
 * SignificanceCommunityCleanUp.cpp
 *
 * Created: 2019-09-11
 * Author: Armin Wiebigke
 */

#include "SignificanceCommunityCleanUp.h"
#include "MergeCommunities.h"

namespace NetworKit {

SignificanceCommunityCleanUp::SignificanceCommunityCleanUp(const Graph &graph,
                                                           std::vector<std::vector<node>> &communities,
                                                           StochasticDistribution &distribution,
                                                           double significanceThreshold,
                                                           double scoreThreshold,
                                                           double minOverlapRatio,
                                                           bool mergeDiscarded)
		: graph(graph),
		  communities(communities),
		  significanceThreshold(significanceThreshold),
		  scoreThreshold(scoreThreshold),
		  minOverlapRatio(minOverlapRatio),
		  mergeDiscarded(mergeDiscarded),
		  maxCommunitySize(0),
		  stochasticDistribution(distribution) {
}

void SignificanceCommunityCleanUp::run() {
	hasRun = false;
	cleanAllCommunities();
	if (mergeDiscarded) {
		mergeDiscardedCommunities();
	}
	hasRun = true;
}

std::string SignificanceCommunityCleanUp::toString() const {
	return "SignificanceCommunityCleanUp";
}

bool SignificanceCommunityCleanUp::isParallel() const {
	return true;
}

void SignificanceCommunityCleanUp::cleanAllCommunities() {
	INFO("Clean ", communities.size(), " communities");
	#pragma omp parallel
	{
		SingleCommunityCleanUp singleCommunityCleanup(graph, stochasticDistribution, scoreThreshold, significanceThreshold, minOverlapRatio);
#pragma omp for schedule(dynamic, 1)
		for (omp_index i = 0; i < static_cast<omp_index>(communities.size()); ++i) {
			const std::vector<node>& inputCommunity = communities[i];
			DEBUG("Clean community ", i, "/", communities.size(), " with size ", inputCommunity.size());
			std::vector<node> cleanedCommunity = singleCommunityCleanup.clean(inputCommunity);

			if (!cleanedCommunity.empty()) {
				communities[i] = std::move(cleanedCommunity);
			} else {
				if (mergeDiscarded) {
					#pragma omp critical
					{
						discardedCommunities.emplace_back(std::move(communities[i]));
					}
				}
				communities[i].clear();
			}
		}
	}

	auto new_end = std::remove_if(communities.begin(), communities.end(),
				      [](const std::vector<node>& c) {
					      return c.empty();
				      });
	communities.erase(new_end, communities.end());

	maxCommunitySize = std::max_element(communities.begin(), communities.end(),
					    [](const std::vector<node> &c1, const std::vector<node> &c2) {
						    return c1.size() < c2.size();
					    })->size();
}


void SignificanceCommunityCleanUp::mergeDiscardedCommunities() {
	INFO("Try to merge ", discardedCommunities.size(), " discarded communities");
	MergeCommunities mergeCommunities(graph, std::move(discardedCommunities), stochasticDistribution,
	                                  significanceThreshold, scoreThreshold, minOverlapRatio,
	                                  2 * maxCommunitySize);
	mergeCommunities.run();
	for (const auto &community : mergeCommunities.getCleanedCommunities()) {
		communities.push_back(community);
	}
}

} /* namespace NetworKit */
