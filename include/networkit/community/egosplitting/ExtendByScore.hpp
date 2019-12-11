/*
 * ExtendByScore.h
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#ifndef EXTENDEDGES_H
#define EXTENDEDGES_H

#include <networkit/community/egosplitting/ExtendEgoNetStrategy.hpp>
#include <networkit/community/egosplitting/EgoSplitting.hpp>
#include <networkit/auxiliary/Timings.hpp>
#include <networkit/auxiliary/SparseVector.hpp>

namespace NetworKit {

class ExtendByScore : public ExtendEgoNetStrategy {
public:
	explicit ExtendByScore(EgoNetData &egoNetData, count maxCandidates,
	                       const Graph &egoGraph, node egoNode);

	void run() override;

	bool isParallel() const override;

	std::string toString() const override;

private:
	using NodeAndScore = std::pair<node, double>;
	SparseVector<double> &nodeScores;
	std::string scoreStrategy;
	SignificanceCalculator &significanceCalculator;
	count outgoingStubs;
	count externalStubs;

	template<typename F>
	std::vector<NodeAndScore> calculateScoresImpl(F calculateScore) const;

	std::vector<NodeAndScore> calculateScores() const;

	void searchForCandidates();

	void takeBestCandidates(std::vector<NodeAndScore> &candidatesAndScores);
};

} /* namespace NetworKit */

#endif //EXTENDEDGES_H
