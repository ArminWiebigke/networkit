/*
 * ParallelTimings.h
 *
 * Created on: 2019-11-05
 * Author: Armin Wiebigke
  */

#ifndef NETWORKIT_PARALLELTIMINGS_HPP
#define NETWORKIT_PARALLELTIMINGS_HPP

#include <unordered_map>
#include <vector>
#include <omp.h>

#include <networkit/auxiliary/Timer.hpp>

namespace NetworKit {

/**
 * A class that provides methods to time parts of an algorithm. Inherit from this class to use it.
 */
class ParallelTimings {
public:
	/**
	 * Get timings for the parts of the algorithm.
	 * @return A map that maps the timer name to its value.
	 */
	std::unordered_map<std::string, double> getTimings() const;

	std::string timingsAsString() const;

protected:
	ParallelTimings();

	void addTime(Aux::Timer &timer, const std::string &name) const;

	void addTimings(const std::unordered_map<std::string, double> &ts,
	                const std::string &prefix = "") const;

	bool timingsEmpty() const;

private:
	mutable std::vector<std::unordered_map<std::string, double>> timings;
};

} /* namespace NetworKit */


#endif //NETWORKIT_PARALLELTIMINGS_HPP
