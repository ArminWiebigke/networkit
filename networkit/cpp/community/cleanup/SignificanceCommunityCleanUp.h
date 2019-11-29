/*
 * SignificanceCommunityCleanUp.h
 *
 * Created: 2019-09-11
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_CLEAN_UP_H
#define NETWORKIT_CLEAN_UP_H

#include "../../base/Algorithm.h"
#include "../../graph/Graph.h"
#include "../../structures/Cover.h"
#include "StochasticSignificance.h"
#include "SingleCommunityCleanUp.h"

namespace NetworKit {

/**
 * @ingroup community
 * An algorithm that aims to improve the quality of (overlapping) communities, e.g. to clean up
 * the communities detected by a (overlapping) community detection algorithm.
 * Based on the statistical significance of the communities.
 */
class SignificanceCommunityCleanUp : public Algorithm {
public:
	using Community = SingleCommunityCleanUp::Community;

	/**
	 * Constructor of the algorithm.
	 * @param	graph	input graph
	 * @param	cover	input cover
	 * @param       distribution The stochastic distribution object to use for stochastic calculations
	 */
	SignificanceCommunityCleanUp(const Graph &graph,
	                             const Cover &cover,
	                             StochasticDistribution &distribution,
	                             double significanceThreshold = 0.1,
	                             double scoreThreshold = 0.1,
	                             double minOverlapRatio = 0.5,
	                             bool mergeDiscarded = true);

	void run() override;

	/**
	 * Get the result cover.
	 * @return A cover containing the cleaned communities
	 */
	Cover getCover();

	std::string toString() const override;

	bool isParallel() const override;

private:

	const Graph &graph;
	const Cover &cover;
	Cover cleanedCommunities;
	std::set<Community> discardedCommunities;
	double significanceThreshold;
	double scoreThreshold;
	double minOverlapRatio;
	const bool mergeDiscarded;
	count maxCommunitySize;

	StochasticDistribution &stochasticDistribution;

	void cleanAllCommunities();

	void mergeDiscardedCommunities();

};
} /* namespace NetworKit */

#endif //NETWORKIT_CLEAN_UP_H
