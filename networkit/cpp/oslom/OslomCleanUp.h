/*
 * OslomCleanUp.h
 *
 * Created: 2019-06-06
 * Author: Armin Wiebigke
 */


#ifndef OSLOMCLEANUP_H
#define OSLOMCLEANUP_H

#include <vector>
#include <string>

#include "../base/Algorithm.h"
#include "../graph/Graph.h"
#include "../structures/Cover.h"

namespace NetworKit {

class OslomCleanUp : public Algorithm {

public:

	OslomCleanUp(const Graph &graph, const Cover &cover);

	OslomCleanUp(const Graph &graph, const Cover &cover, std::vector<std::string> args);

	/**
	 * Run the algorithm.
	 */
	void run() override;

	/**
	 * Get the result cover.
	 * @return The cover
	 */
	Cover getCover();

	/**
	 * Get a string representation of the algorithm.
	 *
	 * @return string representation of algorithm and parameters.
	 */
	std::string toString() const override;

	/**
	 * @return True if algorithm can run multi-threaded.
	 */
	bool isParallel() const override;

private:

	const Graph &graph;
	const Cover &cover;
	std::vector<std::string> args;
	Cover resultCover;
};

} /* namespace NetworKit */


#endif /* OSLOMCLEANUP_H */
