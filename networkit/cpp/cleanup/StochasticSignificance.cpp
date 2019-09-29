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
	assert(kIn <= cOut);
//	count openStubs = extStubs + cOut; // TODO: Include cOut or not?
	count openStubs = extStubs;
	dist.setMaxValue(openStubs);

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
//	double bootInterval = (0.5 + bootRandomness) * exactProb; // TODOe: Use this instead of hypergeom.
	double score = rightCum + bootInterval;
//	double score = rightCum + exactProb; // TODO: Use bootInterval or not?
	assert(score <= 1.001);

	score = std::min(score, 1.);
	bootInterval = std::min(bootInterval, 1. - score);
	bootInterval = std::min(bootInterval, score);
	score = std::max(score, 1e-100);
	return {score, bootInterval};
}

double StochasticSignificance::orderStatistic(double sScore, count externalNodes, count pos) const {
	return dist.rightCumulativeBinomial(sScore, externalNodes, pos);
}

} /* namespace NetworKit */