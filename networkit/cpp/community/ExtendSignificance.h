/*
 * ExtendSignificance.h
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#ifndef EXTENDSIGNIFICANCE_H
#define EXTENDSIGNIFICANCE_H

#include <vector>

#include "../graph/Graph.h"
#include "../structures/NodeMapping.h"
#include "../auxiliary/Timings.h"
#include "../structures/AdjacencyArray.h"
#include "../structures/Partition.h"
#include "../base/Algorithm.h"
#include "EgoSplitting.h"
#include "ExtendScore.h"

namespace NetworKit {

class ExtendSignificance : public ExtendScore, public Timings {

public:
	ExtendSignificance(const EgoNetData &egoNetData,
	                   const Partition &basePartition);

	void run() override;

	std::string toString() const override;

	bool isParallel() const override;

private:
	const Partition &basePartition;

	std::vector<std::pair<node, double>>
	calcSignficance(node externalNode, const Graph &coarseGraph, const NodeMapping &coarseMapping,
	                const std::vector<count> &coarseSizes, double orderedStatPosition) const;
};

} /* namespace NetworKit */

#endif //EXTENDSIGNIFICANCE_H
