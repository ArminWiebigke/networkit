/*
 * ExtendEdges.h
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#ifndef EXTENDEDGES_H
#define EXTENDEDGES_H

#include "ExtendEgoNetStrategy.h"
#include "EgoSplitting.h"
#include "../auxiliary/Timings.h"

namespace NetworKit {

class ExtendEdges : public ExtendEgoNetStrategy {
public:
	explicit ExtendEdges(const EgoNetData &egoNetData, count maxCandidates,
	                     const Graph &egoGraph, node egoNode);

	void run() override;

	bool isParallel() const override;

	std::string toString() const override;

private:
	std::vector<double> &nodeScores;
	double normalizeScore(node v, double score) const;
};

} /* namespace NetworKit */

#endif //EXTENDEDGES_H
