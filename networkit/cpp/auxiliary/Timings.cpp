/*
 * Timings.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#include "Timings.h"

namespace NetworKit {

void Timings::addTime(Aux::Timer &timer, const std::string &name) const {
	timer.stop();
	double elapsed = timer.elapsedNanoseconds();
	timings[name] += elapsed;
	timer.start();
}

void Timings::addTimings(const std::unordered_map<std::string, double> &ts,
		const std::string &prefix) {
	for (auto &t : ts)
		timings[prefix + t.first] += t.second;
}

std::unordered_map<std::string, double> Timings::getTimings() {
	return timings;
}

} /* namespace NetworKit */