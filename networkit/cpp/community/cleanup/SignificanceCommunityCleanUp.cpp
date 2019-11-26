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
                                                           double significanceThreshold,
                                                           double scoreThreshold,
                                                           double minOverlapRatio,
                                                           bool mergeDiscarded)
		: graph(graph),
		  cover(cover),
		  singleCommunityCleanup(graph, scoreThreshold, significanceThreshold, minOverlapRatio),
		  mergeDiscarded(mergeDiscarded),
		  maxCommunitySize(0) {
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
	return false;
}

void SignificanceCommunityCleanUp::cleanAllCommunities() {
	auto inputCommunities = cover.getSubsets();
	index communityId = 0;
	INFO("Clean ", inputCommunities.size(), " communities");
	for (auto const &inputCommunity : inputCommunities) {
		DEBUG("Clean community ", communityId++, "/", inputCommunities.size(), " with size ", inputCommunity.size());
		auto cleanedCommunity = cleanCommunity(inputCommunity);
		if (cleanedCommunity.empty() && mergeDiscarded)
			discardedCommunities.insert(inputCommunity);
		else {
			cleanedCommunities.addSubset(cleanedCommunity);
			maxCommunitySize = std::max(maxCommunitySize, cleanedCommunity.size());
		}
	}
}

SignificanceCommunityCleanUp::Community
SignificanceCommunityCleanUp::cleanCommunity(const Community &inputCommunity) {
	return singleCommunityCleanup.clean(inputCommunity);
}

void SignificanceCommunityCleanUp::mergeDiscardedCommunities() {
	INFO("Try to merge ", discardedCommunities.size(), " discarded communities");
	MergeCommunities mergeCommunities(graph, std::move(discardedCommunities), singleCommunityCleanup,
	                                  2 * maxCommunitySize);
	mergeCommunities.run();
	for (const auto &community : mergeCommunities.getCleanedCommunities()) {
		cleanedCommunities.addSubset(community);
	}
}

} /* namespace NetworKit */
