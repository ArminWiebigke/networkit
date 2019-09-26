/*
 * StochasticDistribution.cpp
 *
 * Created: 2019-09-13
 * Author: Armin Wiebigke
 */

#include <cassert>

#include "StochasticDistribution.h"

namespace NetworKit {

//std::vector<double> StochasticDistribution::data = {};

StochasticDistribution::StochasticDistribution(index maxValue) {
	setMaxValue(maxValue);
}

void StochasticDistribution::setMaxValue(count maxValue) {
	if (maxValue < data.size())
		return;
	data.reserve(maxValue + 1);
	if (data.empty())
		data.push_back(0.0);
	size_t old_size = data.size();
	double f = data.back();
	for (index i = old_size; i <= maxValue; ++i) {
		f += std::log(i);
		data.push_back(f);
	}
}

double StochasticDistribution::binomCoeff(count n, count k) const {
	return std::exp(logBinomCoeff(n, k));
}

double StochasticDistribution::binomialDist(double p, count n, count k) const {
	return std::exp(logBinomCoeff(n, k) + k * std::log(p) + (n - k) * std::log(1 - p));
}

// Calculates the change ratio if we increment bottom to bottom+1
// => returns "top choose bottom+1" / "top choose bottom"
inline double binomialCoeffChangeForIncrement(count top, count bottom) {
	return (double) (top - bottom) / (bottom + 1);
}

// Calculates the change ratio if we decrement bottom to bottom-1
// => returns "top choose bottom-1" / "top choose bottom"
inline double binomialCoeffChangeForDecrement(count top, count bottom) {
	return (double) bottom / (top + 1 - bottom);
}

// Calculates the change ratio of a binomial distribution if k is incremented (k -> k+1) or
// decremented (k -> k-1).
class BinomialChangeRatio {
public:
	BinomialChangeRatio(double p, double n) : n(n) {
		successProbabilityChangeIncrementing = p / (1 - p);
		successProbabilityChangeDecrementing = (1 - p) / p;
	}

	double incrementingK(count k) {
		return successProbabilityChangeIncrementing
		       * binomialCoeffChangeForIncrement(n, k);
	}

	double decrementingK(count k) {
		return successProbabilityChangeDecrementing
		       * binomialCoeffChangeForDecrement(n, k);
	}

private:
	count n;
	double successProbabilityChangeIncrementing;
	double successProbabilityChangeDecrementing;
};

double StochasticDistribution::rightCumulativeBinomial(double p, count n, count k) const {
	assert("Error: n < 0" && n >= 0);
	assert("Error: k < 0" && k >= 0);
	assert("Error: k > n" && k <= n);
	assert("Error: p < 0" && p >= 0);
	if (k == 0)
		return 1;
	if (p - 1 > -1e-11)
		return 1;
	// k is smaller than the expected value, so left cumulative is faster
	if (k < n * p)
		return 1 - leftCumulativeBinomial(p, n, k - 1);

	double startBinom = binomialDist(p, n, k);
	if (startBinom <= 1e-40)
		return 0;

	// Same approach as in rightCumulativeHyper
	double curBinom = 1.;
	double sum = curBinom;
	BinomialChangeRatio changeRatio(p, n);
	for (count x = k; x < n; ++x) {
		curBinom *= changeRatio.incrementingK(x);
		sum += curBinom;
		if (curBinom < precision * sum)
			break;
	}
	assert(startBinom * sum <= 1.001);
	return startBinom * sum;
}

double StochasticDistribution::leftCumulativeBinomial(double p, count n, count k) const {
	if (k == n)
		return 1;
	if (p < 1e-11)
		return 1;
	// k is larger than the expected value, so right cumulative is faster
	if (k > n * p)
		return 1 - rightCumulativeBinomial(p, n, k + 1);

	double startBinom = binomialDist(p, n, k);
	if (startBinom <= 1e-40)
		return 0;

	// Same approach as in rightCumulativeHyper
	double curBinom = 1.;
	double sum = curBinom;
	BinomialChangeRatio changeRatio(p, n);
	for (count x = k; x > 0; --x) {
		curBinom *= changeRatio.decrementingK(x);
		sum += curBinom;
		if (curBinom < precision * sum)
			break;
	}
	assert(startBinom * sum <= 1.001);
	return startBinom * sum;
}

double StochasticDistribution::hypergeometricDist(count N, count K, count n, count k) const {
	return std::max(0., std::exp(logHyper(N, K, n, k)));
}

double StochasticDistribution::logHyper(count N, count K, count n, count k) const {
	return logBinomCoeff(K, k)
	       + logBinomCoeff(N - K, n - k)
	       - logBinomCoeff(N, n);
}

// Calculates the change ratio of a hypergeometric distribution if k is incremented (k -> k+1) or
// decremented (k -> k-1).
class HypergeomChangeRatio {
public:
	HypergeomChangeRatio(count N, count K, count n) : N(N), K(K), n(n) {};

	double incrementingK(count k) {
		return binomialCoeffChangeForIncrement(K, k) *
		       binomialCoeffChangeForDecrement(N - K, n - k);
	}

	double decrementingK(count k) {
		return binomialCoeffChangeForDecrement(K, k) *
		       binomialCoeffChangeForIncrement(N - K, n - k);
	}

private:
	count N;
	count K;
	count n;
};

double StochasticDistribution::rightCumulativeHyper(count N, count K, count n, count k) const {
	assert("Error: N < 0" && N >= 0);
	assert("Error: K < 0" && K >= 0);
	assert("Error: K > N" && K <= N);
	assert("Error: n < 0" && n >= 0);
	assert("Error: k < 0" && k >= 0);
	assert(data.size() > N);
	if (k == 0)
		return 1;
	if (k > n || k > K)
		return 0;
	// If k is smaller than the expected value, calculating the left cumulative is faster
	double expectedValue = (double) K / N * n;
	if (k < expectedValue)
		return (1. - leftCumulativeHyper(N, K, n, k - 1));

	double startProb = hypergeometricDist(N, K, n, k);
	if (startProb <= 1e-40)
		return 0;

	// Sum the probabilities until the sum only changes minimally.
	// To increase the performance, we calculate the change in the probability of the hypergeometric
	// distribution for increasing x instead of calculating each probability directly.
	// At the end, the sum is normalized with the starting probability.
	double curProb = 1.;
	double sum = curProb;
	HypergeomChangeRatio changeRatio(N, K, n);
	for (count x = k; x <= n; ++x) {
		curProb *= changeRatio.incrementingK(x);
		sum += curProb;
		if (curProb < precision * sum)
			break;
	}

	assert(startProb * sum <= 1.001);
	return startProb * sum;
}

double StochasticDistribution::leftCumulativeHyper(count N, count K, count n, count k) const {
	if (N - K - n + k < 0)
		return 0;
	if (k < 0)
		return 0;
	if (k == n)
		return 1;
	// If k is larger than the expected value, calculating the right cumulative is faster
	double expectedValue = (double) K / N * n;
	if (k > expectedValue)
		return (1. - rightCumulativeHyper(N, K, n, k + 1));

	double startProb = hypergeometricDist(N, K, n, k);
	if (startProb <= 1e-40)
		return 0;

	// Same approach as in rightCumulativeHyper
	double curProb = 1.;
	double sum = curProb;
	HypergeomChangeRatio changeRatio(N, K, n);
	for (count x = k; x > 0; --x) {
		curProb *= changeRatio.decrementingK(x);
		sum += curProb;
		if (curProb < precision * sum)
			break;
	}

	assert(startProb * sum <= 1.001);
	return startProb * sum;
}

// Calculates the change ratio of the statistical significance distribution if k is
// incremented (k -> k+1) or decremented (k -> k-1).
class SignificanceChangeRatio {
public:
	SignificanceChangeRatio(count kTotal, count cOut, count extStubs) : kTotal(kTotal), cOut(cOut) {
		count M = extStubs - kTotal;
		MInEdgesConstPart = (M - cOut - kTotal) / 2; // MIn = M - cOut - kTotal + 2*kIn
	}

	double incrementingK(count kIn) {
		return 0.5                    // 2^-kIn
		       * (kTotal - kIn)       // 1/kOut! (-1 in Factorial)
		       / (kIn + 1)            // 1/kIn!  (+1 in Factorial)
		       * (cOut - kIn)         // 1/(cOut - kIn)!  (-1 in Factorial)
		       / (MInEdgesConstPart + (kIn + 1))     // 1/(MIn/2)!  (+1 in Factorial)
				;
	}

	double decrementingK(count kIn) {
		return 2.                           // 2^-kIn
		       / (kTotal - (kIn - 1))       // 1/kOut! (+1 in Factorial)
		       * kIn                        // 1/kIn!  (-1 in Factorial)
		       / (cOut - (kIn - 1))         // 1/(cOut - kIn)!  (+1 in Factorial)
		       * (MInEdgesConstPart + kIn)  // 1/(MIn/2)!  (-1 in Factorial)
				;
	}

private:
	count kTotal;
	count cOut;
	count MInEdgesConstPart;
};

std::pair<double, double>
StochasticDistribution::rightCumulativeStochastic(count kTotal, count kIn, count cOut,
                                                  count extStubs) const {
	if (kIn > cOut || kIn > kTotal)
		return {0.0, 0.0};

//	// Get mode == highest probability
//	// TODO: Can we do this better?
//	count mode = (cOut / double(M + kTotal) * kTotal);
//	mode = std::min(mode, cOut); // this mode is underestimated anyway
//	// TODO: Use mode for speedup?

	// Right cumulative
	double kInProbability = 1.0; // p(kIn) (not normalized)
	double rightCumulativeSum = kInProbability;
	double currentProbability = kInProbability;
	count maxKIn = std::min(kTotal, cOut);
	SignificanceChangeRatio changeRatio(kTotal, cOut, extStubs);
	for (count x = kIn; x < maxKIn; ++x) {
		double ratio = changeRatio.incrementingK(x);
		currentProbability *= ratio;
		rightCumulativeSum += currentProbability;
		assert(ratio != 0.0);
		assert(currentProbability < 1e200);

		double sumChange = currentProbability / rightCumulativeSum;
		if (sumChange < precision)
			break;
	}

	// Left cumulative
	double probabilitySum = rightCumulativeSum;
	currentProbability = kInProbability;
	for (count x = kIn; x > 0; --x) {
		currentProbability *= changeRatio.decrementingK(x);
		probabilitySum += currentProbability;
		assert(x < 1e6);

		double sumChange = currentProbability / probabilitySum;
		if (sumChange < precision)
			break;
	}

	double normalizedKInProbability = kInProbability / probabilitySum;
	double normalizedRightCumulative = rightCumulativeSum / probabilitySum;
	return {normalizedKInProbability, normalizedRightCumulative};
}

} /* namespace NetworKit */