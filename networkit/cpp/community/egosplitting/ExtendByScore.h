/*
 * ExtendByScore.h
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#ifndef EXTENDEDGES_H
#define EXTENDEDGES_H

#include "ExtendEgoNetStrategy.h"
#include "EgoSplitting.h"
#include "../../auxiliary/Timings.h"
#include "../../structures/SparseVector.h"

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
	StochasticSignificance &significance;
	count outgoingStubs;
	count externalStubs;

	double calculateScore(node v, count numEdges) const;

	std::vector<NodeAndScore> calculateScores(const std::vector<node> &candidates) const;

	std::vector<node> searchForCandidates();

	void takeBestCandidates(std::vector<NodeAndScore> &candidatesAndScores);
};

} /* namespace NetworKit */

#endif //EXTENDEDGES_H
