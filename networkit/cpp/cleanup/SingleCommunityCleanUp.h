/*
 * SingleCommunityCleanup.h
 *
 * Created: 2019-09-26
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_SINGLECOMMUNITYCLEANUP_H
#define NETWORKIT_SINGLECOMMUNITYCLEANUP_H


#include "../graph/Graph.h"
#include "StochasticSignificance.h"

namespace NetworKit {

/**
 * Provides a way to clean up communities to improve their quality.
 */
class SingleCommunityCleanUp {
public:
	using Community = std::set<node>;

	explicit SingleCommunityCleanUp(const Graph &graph,
	                                double scoreThreshold = 0.1,
	                                double significanceThreshold = 0.1,
	                                double minOverlapRatio = 0.5);

	/**
	 * Clean a community using statistical significance. If the community is not considered
	 * significant, an empty community is returned.
	 * @param inputCommunity community to clean
	 * @return cleaned community
	 */
	Community clean(const Community &inputCommunity);

private:
	struct ScoreStruct {
		double sScore;
		double bootInterval;
		node candidate;

		ScoreStruct(double s, double b, node u) : sScore(s), bootInterval(b), candidate(u) {};
	};

	const Graph &graph;
	StochasticSignificance stochastic;
	// threshold to decide if a node is significant
	const double significanceThreshold;
	// threshold to discard candidates because they will most likely not be significant
	const double scoreThreshold;
	// threshold to discard communities if they changed too much
	const double minOverlapRatio;

	Community community;
	std::vector<count> edgesToCommunity;         // kIn[u] == number of neighbors of node u that are in the community
	std::vector<int> isInCommunity; // isInCommunity[u] == is node u in the community
	std::vector<int> isCandidate;   // isCandidate[u] == is node u a candidate
	count outgoingCommunityStubs;   // number of edges that have one endpoint inside the community and one outside
	count totalCommunityStubs;      // the sum of the degrees of the nodes in the community
	count externalNodes;
	count externalStubs;

	std::vector<node> getCandidatesAndSetUpCalculation(bool includeNeighbors);

	void reset(const std::vector<node> &candidates);

	Community calculateSignificantNodes(const Community &inputCommunity, bool includeNeighbors);

	Community findSignificantCandidates(std::vector<ScoreStruct> scores) const;

	std::vector<SingleCommunityCleanUp::ScoreStruct>

	calculateCandidateScores(const std::vector<node> &candidates) const;

	void removeWorstNode(std::vector<ScoreStruct> scores);

	bool smallOverlap(const Community &inputCommunity, const Community &cleanedCommunity) const;
};


} /* namespace NetworKit */


#endif //NETWORKIT_SINGLECOMMUNITYCLEANUP_H
