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
#include "../structures/LowToHighDirectedGraph.h"
#include "../auxiliary/Timings.h"
#include "EgoSplitting.h"

namespace NetworKit {

/**
 * Create a partition of the ego-net of a given node.
 * Optionally, the ego-net is extended before the partitioning.
 */
class EgoNetPartition : public CommunityDetectionAlgorithm, public Timings {
	using PartitionFunction = std::function<Partition(const Graph &)>;

public:
	EgoNetPartition(const EgoNetData &egoNetData, node egoNode, Graph egoGraph,
	                PartitionFunction partitionFunction);

	void run() override;

	bool isParallel() const override;

	std::string toString() const override;

	/**
	 * Get the extended ego-net.
	 * @return extended ego-net graph
	 */
	Graph getExtendedEgoGraph() const;

private:
	const LowToHighDirectedGraph &directedG;
	Graph egoGraph;
	NodeMapping &egoMapping;
	node egoNode;
	const PartitionFunction partitionFunction;
	const std::unordered_map<std::string, std::string> &parameters;
	const Cover &groundTruth;
	const EgoNetData &egoNetData;

	void partitionEgoNet();

	void extendEgoNet(const std::string &extendStrategy);

	void removeLowDegreeNodes(count minDegree, count directNeighborsBound);

	Partition createGroundTruthPartition() const;

	count extendIterationsCount() const;

	void extendAndPartition();
};

} /* namespace NetworKit */

#endif /* EGONETPARTITION_H */
