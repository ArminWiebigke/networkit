/*
 * SignificanceCommunityCleanUp.h
 *
 * Created: 2019-09-11
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_CLEAN_UP_H
#define NETWORKIT_CLEAN_UP_H

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Cover.hpp>
#include <networkit/community/cleanup/SignificanceCalculator.hpp>
#include <networkit/community/cleanup/SingleCommunityCleanUp.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * An algorithm that aims to improve the quality of (overlapping) communities, e.g. to clean up
 * the communities detected by a (overlapping) community detection algorithm.
 * Based on the statistical significance of the communities.
 */
class SignificanceCommunityCleanUp : public Algorithm {
public:
	/**
	 * Constructor of the algorithm.
	 * @param	graph	input graph
	 * @param	cover	input cover
	 * @param       distribution The stochastic distribution object to use for stochastic calculations
	 */
	SignificanceCommunityCleanUp(const Graph &graph,
	                             std::vector<std::vector<node>> &communities,
	                             StochasticDistribution &distribution,
	                             double significanceThreshold = 0.1,
	                             double scoreThreshold = 0.1,
	                             double minOverlapRatio = 0.5,
	                             bool mergeDiscarded = true);

	void run() override;

	std::string toString() const override;

	bool isParallel() const override;

private:

	const Graph &graph;
	std::vector<std::vector<node>> &communities;
	std::vector<std::vector<node>> discardedCommunities;
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
