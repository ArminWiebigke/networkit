/*
 * SignificanceCommunityCleanUp.h
 *
 * Created: 2019-09-11
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_CLEAN_UP_H
#define NETWORKIT_CLEAN_UP_H

#include "../base/Algorithm.h"
#include "../graph/Graph.h"
#include "../structures/Cover.h"
#include "StochasticSignificance.h"

namespace NetworKit {


class SignificanceCommunityCleanUp : public Algorithm {
	using Community = std::set<index>;

public:

	SignificanceCommunityCleanUp(const Graph &graph, const Cover &cover);

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
	// threshold to discard communities if they changed too much
	double gamma;
	std::vector<count> kIn;
	std::vector<bool> inComm;
	StochasticSignificance stochastic;
	double sigThreshold = 0.1;

	struct ScoreStruct {
		double sScore;
		double bootInterval;
		node u;

		ScoreStruct(double s, double b, node u) : sScore(s), bootInterval(b), u(u) {};
	};

	std::vector<Community> discardedComms;

	Community cleanCommunity(const Community &community);

	void checkComms();

	void mergeDiscarded();

	SignificanceCommunityCleanUp::Community checkCommSignificance(Community community, bool checkNeighbors);

	Community significantCandidates(const std::vector<ScoreStruct>& scores, count externalNodes);

};
} /* namespace NetworKit */

#endif //NETWORKIT_CLEAN_UP_H
