/*
 * MergeCommunities.h
 *
 * Created: 2019-09-26
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_MERGECOMMUNITIES_H
#define NETWORKIT_MERGECOMMUNITIES_H

#include "../../graph/Graph.h"
#include "SingleCommunityCleanUp.h"
#include "../../structures/Partition.h"
#include "../../base/Algorithm.h"

namespace NetworKit {

/**
 * Merge not significant communities to find significant communities, using the statistical
 * significance.
 */
class MergeCommunities : public Algorithm {
public:
	using Community = SingleCommunityCleanUp::Community;

	MergeCommunities(const Graph &graph,
	                 std::set<Community> discardedCommunities,
	                 StochasticDistribution &stochasticDistribution,
	                 double significanceThreshold = 0.1,
	                 double scoreThreshold = 0.1,
	                 double minOverlapRatio = 0.5,
	                 count maxCommunitySize = none);

	void run() override;

	std::vector<Community> getCleanedCommunities();

	std::string toString() const override;

	bool isParallel() const override;

private:
	const Graph &graph;
	std::set<Community> discardedCommunities;
	StochasticDistribution &stochasticDistribution;
	SignificanceCalculator significanceCalculator;
	double significanceThreshold;
	double scoreThreshold;
	double minOverlapRatio;
	std::vector<Community> cleanedCommunities;
	Graph discardedCommunitiesGraph;
	Partition mergedCommunities;
	std::vector<std::set<node>> coarseToFineMapping;
	std::vector<count> outgoingGroupStubs;
	std::vector<count> totalGroupStubs;
	count totalStubs;
	const count maxCommunitySize;

	void createDiscardedCommunitiesGraph();

	void tryToMergeCommunities();

	void checkMergedCommunities();

	bool tryLocalMove(node u);
};


} /* namespace NetworKit */


#endif //NETWORKIT_MERGECOMMUNITIES_H
