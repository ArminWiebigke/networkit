/*
 * StochasticSignificance.h
 *
 * Created: 2019-09-13
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_STOCHASTICSIGNIFICANCE_H
#define NETWORKIT_STOCHASTICSIGNIFICANCE_H

#include "StochasticDistribution.h"
#include "../Globals.h"

namespace NetworKit {

/**
 * This class calculates the statistical significance of a node to a community.
 */
class StochasticSignificance {
public:

	explicit StochasticSignificance(count maxValue);

	/**
	 * Calculate the s-score
	 * @param kIn Number of edges between node and community
	 * @param cOut Number of outgoing stubs from the community
	 * @param extStubs Number of stubs in the rest of the graph (without the node and the community)
	 * @param k Degree of the node
	 * @return a pair (s-score, boot interval)
	 */
	double rScore(count k, count kIn, count cOut, count extStubs) const;

	/**
	 * Calculate the order statistic
	 * @param rScore the s-Score of the candidate
	 * @param externalNodes the number of external nodes
	 * @param pos the position of the candidate
	 * @return
	 */
	double orderStatistic(double rScore, count externalNodes, count pos) const;

private:
	mutable StochasticDistribution dist;
};

} /* namespace NetworKit */

#endif //NETWORKIT_STOCHASTICSIGNIFICANCE_H
