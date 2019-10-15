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
                                                           double minOverlapRatio)
		: graph(graph),
		  cover(cover),
		  singleCommunityCleanup(graph, scoreThreshold, significanceThreshold, minOverlapRatio) {
}

void SignificanceCommunityCleanUp::run() {
	hasRun = false;
	cleanedCommunities = Cover(graph.upperNodeIdBound());
	cleanAllCommunities();
	mergeDiscardedCommunities();
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
	for (const Community &inputCommunity : inputCommunities) {
		auto cleanedCommunity = cleanCommunity(inputCommunity);
		if (cleanedCommunity.empty())
			discardedCommunities.insert(inputCommunity);
		else
			cleanedCommunities.addSubset(cleanedCommunity);
	}
}

SignificanceCommunityCleanUp::Community
SignificanceCommunityCleanUp::cleanCommunity(const Community &inputCommunity) {
	return singleCommunityCleanup.clean(inputCommunity);
}

void SignificanceCommunityCleanUp::mergeDiscardedCommunities() {
	MergeCommunities mergeCommunities(graph, std::move(discardedCommunities), singleCommunityCleanup);
	mergeCommunities.run();
	for (const auto &community : mergeCommunities.getCleanedCommunities()) {
		cleanedCommunities.addSubset(community);
	}
}

} /* namespace NetworKit */
