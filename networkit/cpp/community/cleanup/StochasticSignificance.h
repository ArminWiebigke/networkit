/*
 * StochasticSignificance.h
 *
 * Created: 2019-09-13
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_STOCHASTICSIGNIFICANCE_H
#define NETWORKIT_STOCHASTICSIGNIFICANCE_H

#include <random>

#include "StochasticDistribution.h"
#include "../../Globals.h"

namespace NetworKit {

/**
 * This class calculates the statistical significance of a node to a community.
 */
class StochasticSignificance {
public:

	/**
	 * Constructor
	 * @param maxValue maximum value that can be used as an argument for the calculations
	 */
	explicit StochasticSignificance(const StochasticDistribution& dist);

	/**
	 * Calculate the r-score
	 * @param kIn Number of edges between node and community
	 * @param cOut Number of outgoing stubs from the community
	 * @param extStubs Number of stubs in the rest of the graph (without the node and the community)
	 * @param k Degree of the node
	 * @return a pair (s-score, boot interval)
	 */
	double rScore(count k, count kIn, count cOut, count extStubs);

	/**
	 * Calculate the order statistic (s-score)
	 * @param rScore the r-score of the candidate
	 * @param externalNodes the number of external nodes
	 * @param pos the position of the candidate
	 * @return
	 */
	double orderStatistic(double rScore, count externalNodes, count pos);

private:
	const StochasticDistribution& dist;
	std::mt19937_64 rng;
	std::uniform_real_distribution<double> random_distribution;

	void ensureMaxValue(count maxValue) {
		if (dist.maxValue() < maxValue)
			throw std::runtime_error("Maximum value of the distribution is not high enough.");
	}
};

} /* namespace NetworKit */

#endif //NETWORKIT_STOCHASTICSIGNIFICANCE_H
