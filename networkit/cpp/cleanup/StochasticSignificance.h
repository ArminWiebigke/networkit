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
	std::pair<double, double> sScore(count k, count kIn, count cOut, count extStubs);

	/**
	 * Calculate the order statistic
	 * @param sScore the s-Score of the candidate
	 * @param externalNodes the number of external nodes
	 * @param pos the position of the candidate
	 * @return
	 */
	double orderStatistic(double sScore, count externalNodes, count pos);

private:
	StochasticDistribution dist;
};

} /* namespace NetworKit */

#endif //NETWORKIT_STOCHASTICSIGNIFICANCE_H
