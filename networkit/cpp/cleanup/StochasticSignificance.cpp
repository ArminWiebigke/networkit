/*
 * StochasticSignificance.cpp
 *
 * Created: 2019-09-13
 * Author: Armin Wiebigke
 */

#include <algorithm>

#include "StochasticSignificance.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

StochasticSignificance::StochasticSignificance(count maxValue) : dist(maxValue) {

}

std::pair<double, double>
StochasticSignificance::sScore(count k, count kIn, count cOut, count extStubs) const {
	// TODO: Include k in extStubs? currently is not
//	extStubs -= k;
	assert(kIn <= cOut);
	count openStubs = extStubs + cOut;

	double exactProb = 0, rightCum = 0;
	if (k * k < openStubs) {
		// Use approximation with hypergeometric distribution
		exactProb = dist.hypergeometricDist(openStubs, cOut, k, kIn);
		rightCum = dist.rightCumulativeHyper(openStubs, cOut, k, kIn + 1);
	} else {
		// Calculate the probability using the original distribution
		// extInternal = extStubs - (cOut - (k - kIn) + kIn)
		// extInternal = extStubs - cOut + k - 2*kIn
		// exactProb = A * 2^-kIn / ((k-kIn)! * kIn! * (cOut-kIn)! * extStubs!
		// kIn += 1 => /2, *(k-kIn), /(kIn+1), *(cOut-kIn), *(extStubs-1)*(extStubs-2)
		// TODO: Implement caluclation
//		throw std::runtime_error("Exact calculation not implemented!");
		exactProb = dist.hypergeometricDist(openStubs, cOut, k, kIn);
		rightCum = dist.rightCumulativeHyper(openStubs, cOut, k, kIn + 1);
	}
//	double score = rightCum + bootInterval;
	double score = rightCum + exactProb;
	assert(score <= 1.001);
	score = std::min(score, 1.);
	score = std::max(score, 1e-100);

	double bootRandomness = Aux::Random::real(-0.5, 0.5) * 1e-6;
	double bootInterval = (0.5 + bootRandomness) * exactProb;
	bootInterval = std::min(bootInterval, 1. - score);
	bootInterval = std::min(bootInterval, score);
	return {score, bootInterval};
}

double StochasticSignificance::orderStatistic(double sScore, count externalNodes, count pos) const {
	return dist.rightCumulativeBinomial(sScore, externalNodes, pos);
}


} /* namespace NetworKit */