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
 * Stores the cumulative logarithmus to efficiently calculate binomial and hypergeometic
 * distributions.
 */
class StochasticDistribution {
public:
	explicit StochasticDistribution(index size);

	/**
	 * Calculate the binomial coefficient "n choose k". Returns a floating point number that may
	 * slightly differ from the exact integer result.
	 * @return binomial coefficient
	 */
	double binomCoeff(count n, count k);

	/**
	 * Calculate the binomial distribution for k success.
	 * @param p probability of a success
	 * @param n number of trials
	 * @param k number of successful trials
	 * @return probability of k successes
	 */
	double binomialDist(double p, count n, count k);

	/**
	 * Calculate the cumulative binomial distribution that there are k or more success.
	 * @param p probability of a success
	 * @param n number of trials
	 * @param k number of successful trials
	 * @return probability of k or more successes
	 */
	double rightCumulativeBinomial(double p, count n, count k);

	/**
	 * Calculate the cumulative binomial distribution that there are k or less success.
	 * @param p probability of a success
	 * @param n number of trials
	 * @param k number of successful trials
	 * @return probability of k or less successes
	 */
	double leftCumulativeBinomial(double p, count n, count k);

	/**
	 * Calculate a hypergeometric distribution
	 * @param N number of elements
	 * @param K number of successful elements
	 * @param n number of trials
	 * @param k number of successes
	 * @return probability of k successes
	 */
	double hypergeometricDist(count N, count K, count n, count k);

	/**
	 * Calculate the cumulative hypergeometric distribution that there are k or more success.
	 * @param N number of elements
	 * @param K number of successful elements
	 * @param n number of trials
	 * @param k number of successes
	 * @return probability of k or more successes
	 */
	double rightCumulativeHyper(count N, count K, count n, count k);

	/**
	 * Calculate the cumulative hypergeometric distribution that there are k or less success.
	 * @param N number of elements
	 * @param K number of successful elements
	 * @param n number of trials
	 * @param k number of successes
	 * @return probability of k or less successes
	 */
	double leftCumulativeHyper(count N, count K, count n, count k);

	double oslomDist(count k, count kIn, count cOut, count openStubs);

private:
//	static std::vector<double> data;
	std::vector<double> data;
	static constexpr double precision = 1e-6;

	double logHyper(count N, count K, count n, count k);

	/**
     * @return natural logarithm of the binomial coefficient "n choose k"
     */
	inline double logBinomCoeff(count n, count k) {
		return data[n] - data[n - k] - data[k];
	};
};

} /* namespace NetworKit */

#endif //NETWORKIT_STOCHASTICDISTRIBUTION_H
