/*
 * SignificanceCommunityCleanUp.h
 *
 * Created: 2019-09-11
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_CLEAN_UP_H
#define NETWORKIT_CLEAN_UP_H

#include "../base/Algorithm.h"
#include "../graph/Graph.h"
#include "../structures/Cover.h"
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
	using Community = std::set<index>;

	/**
	 * Constructor of the algorithm.
	 * @param	graph	input graph
	 * @param	cover	input cover
	 */
	SignificanceCommunityCleanUp(const Graph &graph,
	                             const Cover &cover,
	                             double significanceThreshold = 0.1,
	                             double scoreThreshold = 0.1,
	                             double minOverlapRatio = 0.5);

	/**
	 * Run the algorithm.
	 */
	void run() override;

	/**
	 * Get the result cover.
	 * @return cover containing the cleaned communities
	 */
	Cover getCover();

	/**
	 * Get a string representation of the algorithm.
	 *
	 * @return string representation of algorithm and parameters.
	 */
	std::string toString() const override;

	/**
	 * @return True if algorithm can run multi-threaded.
	 */
	bool isParallel() const override;

private:

	const Graph &graph;
	const Cover &cover;
	Cover cleanedCommunities;
	std::set<Community> discardedCommunities;
	SingleCommunityCleanUp singleCommunityCleanup;

	void cleanAllCommunities();

	Community cleanCommunity(const Community &inputCommunity);

	void mergeDiscardedCommunities();

};
} /* namespace NetworKit */

#endif //NETWORKIT_CLEAN_UP_H
