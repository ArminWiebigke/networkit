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



/**
 * @ingroup community
 * An algorithm that aims to improve the quality of (overlapping) communities, e.g. to clean up
 * the communities detected by a (overlapping) community detection algorithm.
 * Based on the statistical significance of the communities.
 */
class SignificanceCommunityCleanUp : public Algorithm {
	using Community = std::set<index>;

public:

	/**
	 * Constructor of the algorithm.
	 * @param	graph	input graph
	 * @param	cover	input cover
	 */
	SignificanceCommunityCleanUp(const Graph &graph, const Cover &cover);

	/**
	 * Run the algorithm.
	 */
	void run() override;

	/**
	 * Get the result cover.
	 * @return cover containing the cleaned communities
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
	const double gamma;
	// threshold to decide if a node is significant
	static constexpr double sigThreshold = 0.1;
	// threshold to discard candidates because the will most likely not be significant
	static constexpr double scoreThreshold = 0.1;
	std::vector<count> kIn;
	std::vector<int> inComm;
	StochasticSignificance stochastic;

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

	void init();
};
} /* namespace NetworKit */

#endif //NETWORKIT_CLEAN_UP_H
