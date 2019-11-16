/*
 * StochasticSignificance.cpp
 *
 * Created: 2019-09-13
 * Author: Armin Wiebigke
 */

#include <algorithm>

#include "StochasticSignificance.h"
#include "../../auxiliary/Random.h"

namespace NetworKit {

StochasticSignificance::StochasticSignificance(count maxValue) : dist(maxValue) {

}

double
StochasticSignificance::rScore(count k, count kIn, count cOut, count extStubs) const {
	assert(kIn <= cOut);
	count openStubs = extStubs + cOut;
	assert(openStubs >= k);
	ensureMaxValue(openStubs);

	double exactProb = 0, rightCum = 0;
	bool lowSelfLoopProbability = k * k < openStubs;
	if (lowSelfLoopProbability) {
		// Use approximation with hypergeometric distribution
		exactProb = dist.hypergeometricDist(openStubs, cOut, k, kIn);
		rightCum = dist.rightCumulativeHyper(openStubs, cOut, k, kIn + 1);
	} else {
		// Calculate the probability using the original distribution
		std::tie(exactProb, rightCum) = dist.rightCumulativeStochastic(k, kIn, cOut, extStubs);
		rightCum -= exactProb;
	}

	double bootRandomness = Aux::Random::real(-0.5, 0.5) * 1e-6;
	double bootInterval = (0.5 + bootRandomness) * dist.hypergeometricDist(openStubs, cOut, k, kIn);
//	double bootInterval = (0.5 + bootRandomness) * exactProb; // TODO: Use this instead of hypergeom.
	double score = rightCum + bootInterval;
//	double score = rightCum + exactProb; // TODO: Use bootInterval or not?
//	assert(score <= 1.001);

	score = std::min(score, 1.);
	score = std::max(score, 1e-100);
	return score;
}

double StochasticSignificance::orderStatistic(double rScore, count externalNodes, count pos) const {
	assert(pos > 0);
	ensureMaxValue(externalNodes);
	return dist.rightCumulativeBinomial(rScore, externalNodes, pos);
}

void StochasticSignificance::ensureMaxValue(count maxValue) const {
	if (dist.maxValue() < maxValue)
		dist.setMaxValue(maxValue);
}

} /* namespace NetworKit */