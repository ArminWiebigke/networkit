/*
 * EgoNetPartition.h
 *
 * Created: 2019-06-17
 * Author: Armin Wiebigke
 */

#ifndef EGONETPARTITION_H
#define EGONETPARTITION_H

#include <unordered_map>
#include <string>

#include "CommunityDetectionAlgorithm.h"
#include "../graph/Graph.h"
#include "../auxiliary/Timer.h"
#include "../structures/NodeMapping.h"
#include "../structures/Cover.h"
#include "../structures/AdjacencyArray.h"
#include "../auxiliary/Timings.h"
#include "EgoSplitting.h"

namespace NetworKit {

class EgoNetPartition : public CommunityDetectionAlgorithm, public Timings {
	using PartitionFunction = std::function<Partition(const Graph &)>;

public:
	EgoNetPartition(const EgoNetData &egoNetData,
	                const PartitionFunction &partitionFunction);

	bool isParallel() const override;

	void run() override;

	/**
	 * Get a string representation of the algorithm.
	 *
	 * @return string representation of algorithm and parameters.
	 */
	std::string toString() const override;

private:
	const AdjacencyArray &directedG;
	Graph &egoGraph;
	NodeMapping &egoMapping;
	node egoNode;
	const PartitionFunction &partitionFunction;
	const std::unordered_map<std::string, std::string> &parameters;
	const Cover &groundTruth;
	const EgoNetData &egoNetData;
	int it_char;

	void partitionEgoNet();

	void extendEgoNet(const std::string &extendStrategy);

	void removeLowDegreeNodes(count minDegree, count directNeighborsCnt) const;

	Partition createGroundTruthPartition() const;

};

} /* namespace NetworKit */


#endif /* EGONETPARTITION_H */
