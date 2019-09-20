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

StochasticDistribution::StochasticDistribution(index size) {
//	if (size < data.size())
//		return;
//	data.reserve(size + 1);
//	if (data.empty())
//		data.push_back(0.0);
//	size_t old_size = data.size();
//	double f = data.back();
//	for (index i = old_size; i <= size; ++i) {
//		f += std::log(i);
//		data.push_back(f);
//	}
	data.push_back(0.0);
	double f = 0.0;
	for (index i = 1; i <= size; ++i) {
		f += std::log(i);
		data.push_back(f);
	}
}

double StochasticDistribution::binomCoeff(count n, count k) {
	return std::exp(logBinomCoeff(n, k));
}

double StochasticDistribution::binomialDist(double p, count n, count k) {
	return std::exp(logBinomCoeff(n, k) + k * std::log(p) + (n - k) * std::log(1 - p));
}

double StochasticDistribution::rightCumulativeBinomial(double p, count n, count k) {
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
	double probChange = p / (1 - p);
	for (count x = k + 1; x <= n; ++x) {
		// TODO: Precision problems?
		curBinom *= double(n + 1 - x) / x; // binomial coefficient change
		curBinom *= probChange; // success probability change
		sum += curBinom;
		if (curBinom < precision * sum)
			break;
	}
	return startBinom * sum;
}

double StochasticDistribution::leftCumulativeBinomial(double p, count n, count k) {
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
	double probChange = (1 - p) / p;
	for (count x = k; x > 0; --x) {
		curBinom *= x / double(n + 1 - x); // binomial coefficient change
		curBinom *= probChange; // success probability change
		sum += curBinom;
		if (curBinom < precision * sum)
			break;
	}
	return startBinom * sum;
}

double StochasticDistribution::hypergeometricDist(count N, count K, count n, count k) {
	double logH = logHyper(N, K, n, k);
	return std::max(0., std::exp(logH));
}

double StochasticDistribution::logHyper(count N, count K, count n, count k) {
	double a = logBinomCoeff(K, k);
	double b = logBinomCoeff(N - K, n - k);
	double c = logBinomCoeff(N, n);
	return a
	       + b
	       - c;
}

double StochasticDistribution::rightCumulativeHyper(count N, count K, count n, count k) {
	assert("Error: N < 0" && N >= 0);
	assert("Error: K < 0" && K >= 0);
	assert("Error: K > N" && K <= N);
	assert("Error: n < 0" && n >= 0);
	assert("Error: k < 0" && k >= 0);
	if (k == 0)
		return 1;
	if (k > n || k > K)
		return 0;
	// k is smaller than the expected value, so left cumulative is faster
	if (k < double(K) / N * n)
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
	for (count x = k + 1; x <= n; ++x) {
		curProb *= double(K + 1 - x) / x; // (K choose x-1) -> (K choose x)
		count oldChoose = n - (x - 1);
		curProb *= double(oldChoose) / (N - K + 1 - oldChoose); // (N-K choose n-(x-1)) -> (N-k choose n-x)
		sum += curProb;
		if (curProb < precision * sum)
			break;
	}

	return startProb * sum;
}

double StochasticDistribution::leftCumulativeHyper(count N, count K, count n, count k) {
	if (N - K - n + k < 0)
		return 0;
	if (k < 0)
		return 0;
	if (k == n)
		return 1;
	// k is larger than the expected value, so right cumulative is faster
	if (k > double(K) / N * n)
		return (1. - rightCumulativeHyper(N, K, n, k + 1));

	double startProb = hypergeometricDist(N, K, n, k);
	if (startProb <= 1e-40)
		return 0;

	// Same approach as in rightCumulativeHyper
	double curProb = 1.;
	double sum = curProb;
	for (count x = k; x > 0; --x) {
		curProb *= x / double(K + 1 - x); // (K choose x) -> (K choose x-1)
		count newChoose = n - (x - 1);
		curProb *= double(N - K + 1 - newChoose) / (newChoose); // (N-K choose n-x) -> (N-k choose n-(x-1))
		sum += curProb;
		if (curProb < precision * sum)
			break;
	}

	return startProb * sum;
}


} /* namespace NetworKit */