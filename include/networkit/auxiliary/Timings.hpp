/*
 * Timings.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#ifndef TIMINGS_H
#define TIMINGS_H

#include <unordered_map>
#include <string>

#include <networkit/auxiliary/Timer.hpp>

namespace NetworKit {

/**
 * A class that provides methods to time parts of an algorithm.
 */
class Timings {
public:

	/**
	 * Get timings for the parts of the algorithm.
	 * @return A map that maps the timer name to its value.
	 */
	std::unordered_map<std::string, double> getTimings();

protected:
	void addTime(Aux::Timer &timer, const std::string &name) const;

	void addTimings(const std::unordered_map<std::string, double> &ts,
			const std::string &prefix = "");

	mutable std::unordered_map<std::string, double> timings;
};

} /* namespace NetworKit */

#endif //TIMINGS_H
