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
		node nodeId;

		explicit ScoreStruct(double score, double boot, node nodeId)
				: sScore(score), bootInterval(boot), nodeId(nodeId) {};
	};

	const Graph &graph;
	StochasticSignificance stochastic;
	// threshold to decide if a node is significant
	const double significanceThreshold;
	// threshold to discard candidates because they will most likely not be significant
	const double scoreThreshold;
	// threshold to discard communities if they changed too much
	const double minOverlapRatio;

	Community originalCommunity;
	Community community;
	std::vector<node> candidates;
	std::vector<count> edgesToCommunity;
	std::vector<int> isInCommunity;
	std::vector<int> isInOriginalCommunity;
	std::vector<int> isCandidate;
	count outgoingCommunityStubs;   // number of edges that have one endpoint inside the community and one outside
	count totalCommunityStubs;      // the sum of the degrees of the nodes in the community
	count externalNodes;
	count externalStubs;

	void getCandidatesAndSetUpCalculation(bool onlyUseOriginalCommunity);

	void reset();

	Community
	calculateSignificantNodes(const Community &inputCommunity, bool onlyUseOriginalCommunity);

	std::vector<node> findSignificantCandidates(std::vector<ScoreStruct> scores) const;

	std::vector<SingleCommunityCleanUp::ScoreStruct>

	calculateCandidateScores() const;

	void removeWorstNode(std::vector<ScoreStruct> internalScores);

	bool smallOverlap(const Community &inputCommunity, const Community &cleanedCommunity) const;

	std::vector<ScoreStruct> calculateInternalScores();
};

} /* namespace NetworKit */

#endif //NETWORKIT_SINGLECOMMUNITYCLEANUP_H
