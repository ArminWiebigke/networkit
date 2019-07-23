/*
 * ExtendEdges.h
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#ifndef EXTENDEDGES_H
#define EXTENDEDGES_H

#include "ExtendScore.h"
#include "EgoSplitting.h"
#include "../auxiliary/Timings.h"

namespace NetworKit {

class ExtendEdges : public ExtendScore {
public:
	explicit ExtendEdges(const EgoNetData &egoNetData, count maxCandidates);

	void run() override;

	bool isParallel() const override;

	std::string toString() const override;

private:
	std::vector<double> &nodeScores;
	double normalizeScore(node v, double score) const;
};

} /* namespace NetworKit */

#endif //EXTENDEDGES_H
