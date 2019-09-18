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
StochasticSignificance::sScore(count kIn, count cOut, count extStubs, count k) {
	// TODO: Include k in extStubs? currently is not
//	extStubs -= k;
	count openStubs = extStubs + cOut;
	double rand = 1e-6 * (Aux::Random::real(-0.5, 0.5));
	double hyper = dist.hypergeometricDist(openStubs, cOut, k, kIn);
	double bootInterval = (0.5 + rand) * hyper;
	// TODO: Is bootInterval correct if k*k > openStubs
	double rightCum;
	if (kIn == k || kIn == cOut) {
		rightCum = 0;
	} else if (k * k < openStubs) {
		rightCum = dist.rightCumulativeHyper(openStubs, cOut, k, kIn + 1);
		// TODO: What happens if k > cOut
	} else {
		rightCum = dist.rightCumulativeHyper(openStubs, cOut, k, kIn + 1);
//		throw std::runtime_error("Exact calculation not implemented!");
	}
//	double score = rightCum + bootInterval;
	double score = rightCum + hyper;

	if (score > 1)
		score = 1;
	bootInterval = std::min(bootInterval, 1. - score);
	bootInterval = std::min(bootInterval, score);
	return {std::max(score, 1e-100), bootInterval};
}

double StochasticSignificance::orderStatistic(double sScore, count externalNodes, count pos) {
	return dist.rightCumulativeBinomial(sScore, externalNodes, pos);
}


} /* namespace NetworKit */