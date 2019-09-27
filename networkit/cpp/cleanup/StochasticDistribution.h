/*
 * StochasticDistribution.h
 *
 * Created: 2019-09-13
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_STOCHASTICDISTRIBUTION_H
#define NETWORKIT_STOCHASTICDISTRIBUTION_H

#include <vector>

#include "../Globals.h"

namespace NetworKit {

/**
 * Stores the cumulative logarithmus to efficiently calculate stochastic distributions.
 */
class StochasticDistribution {
public:
	explicit StochasticDistribution(index maxValue);

	void setMaxValue(count maxValue);

	/**
	 * Calculate the binomial coefficient "n choose k". Returns a floating point number that may
	 * slightly differ from the exact integer result.
	 * @return binomial coefficient
	 */
	double binomCoeff(count n, count k) const;

	/**
	 * Calculate the binomial distribution for k success.
	 * @param p probability of a success
	 * @param n number of trials
	 * @param k number of successful trials
	 * @return probability of k successes
	 */
	double binomialDist(double p, count n, count k) const;

	/**
	 * Calculate the cumulative binomial distribution that there are k or more success.
	 * @param p probability of a success
	 * @param n number of trials
	 * @param k number of successful trials
	 * @return probability of k or more successes
	 */
	double rightCumulativeBinomial(double p, count n, count k) const;

	/**
	 * Calculate the cumulative binomial distribution that there are k or less success.
	 * @param p probability of a success
	 * @param n number of trials
	 * @param k number of successful trials
	 * @return probability of k or less successes
	 */
	double leftCumulativeBinomial(double p, count n, count k) const;

	/**
	 * Calculate a hypergeometric distribution
	 * @param N number of elements
	 * @param K number of successful elements
	 * @param n number of trials
	 * @param k number of successes
	 * @return probability of k successes
	 */
	double hypergeometricDist(count N, count K, count n, count k) const;

	/**
	 * Calculate the cumulative hypergeometric distribution that there are k or more success.
	 * @param N number of elements
	 * @param K number of successful elements
	 * @param n number of trials
	 * @param k number of successes
	 * @return probability of k or more successes
	 */
	double rightCumulativeHyper(count N, count K, count n, count k) const;

	/**
	 * Calculate the cumulative hypergeometric distribution that there are k or less success.
	 * @param N number of elements
	 * @param K number of successful elements
	 * @param n number of trials
	 * @param k number of successes
	 * @return probability of k or less successes
	 */
	double leftCumulativeHyper(count N, count K, count n, count k) const;

	/**
	 * Calculate the probability that a node has kIn or more edges to the community, according to
	 * the null model
	 * @param kTotal degree of the node
	 * @param kIn number of edges between node and community
	 * @param cOut number of outgoing stubs from the community
	 * @param extStubs number of stubs in the rest of the graph
	 * @return a pair (p(x = kIn), p(x >= kIn))
	 */
	std::pair<double, double> rightCumulativeStochastic(count kTotal, count kIn, count cOut,
	                                                    count extStubs) const;

private:
//	static std::vector<double> data;
	std::vector<double> data;
	static constexpr double precision = 1e-6;

	double logHyper(count N, count K, count n, count k) const;

	/**
     * @return natural logarithm of the binomial coefficient "n choose k"
     */
	inline double logBinomCoeff(count n, count k) const {
		return data[n] - data[n - k] - data[k];
	};

};

} /* namespace NetworKit */

#endif //NETWORKIT_STOCHASTICDISTRIBUTION_H
