/*
 * EgoSplitting.h
 *
 * Created: 2018-12-11
 * Author: Armin Wiebigke
 */

#ifndef EGOSPLITTING_H
#define EGOSPLITTING_H

#include <unordered_map>
#include <functional>
#include <ostream>

#include "../Globals.h"
#include "../base/Algorithm.h"
#include "../structures/Cover.h"
#include "../structures/AdjacencyArray.h"
#include "../structures/NodeMapping.h"
#include "../auxiliary/Timer.h"
#include "../auxiliary/Timings.h"

namespace NetworKit {

/**
 * Ego Splitting is a framework to detect overlapping communities.
 * https://dl.acm.org/citation.cfm?id=3098054
 */
class EgoSplitting : public Algorithm, public Timings {
	using PartitionFunction = std::function<Partition(const Graph &)>;

public:
	/**
	 * Construct an instance of this algorithm for the input graph.
	 *
	 * @param[in]	G   input graph
	 * @param[in]   localClusterAlgo    algorithm to cluster the ego-net
	 * @param[in]   globalClusterAlgo   algorithm to cluster the persona graph
	 */
	explicit EgoSplitting(const Graph &G);

	/**
	 * Construct an instance of this algorithm for the input graph.
	 *
	 * @param[in]	G   input graph
	 * @param[in]   localClusterAlgo    algorithm to cluster the ego-net
	 * @param[in]   globalClusterAlgo   algorithm to cluster the persona graph
	 */
	EgoSplitting(const Graph &G,
	             PartitionFunction clusterAlgo);

	/**
	 * Construct an instance of this algorithm for the input graph.
	 *
	 * @param[in]	G   input graph
	 * @param[in]   localClusterAlgo    algorithm to cluster the ego-net
	 * @param[in]   globalClusterAlgo   algorithm to cluster the persona graph
	 */
	EgoSplitting(const Graph &G,
	             PartitionFunction localClusterAlgo,
	             PartitionFunction globalClusterAlgo);

	/**
	 * Detect communities.
	 */
	void run() override;

	/**
	 * Returns the result of the run method or throws an error, if the algorithm hasn't run yet.
	 * @return partition of the node set
	 */
	Cover getCover();

	/**
	 * Get a string representation of the algorithm.
	 *
	 * @return string representation of algorithm and parameters.
	 */
	std::string toString() const override;

	/**
	 * Get additional information about the execution.
	 * @return A map that maps a name to a value.
	 */
	std::unordered_map<std::string, double> getExecutionInfo();

	/**
	 * Get the partitions of the ego-nets. A partition maps a node to its partition ID.
	 * @return A vector of the partitions.
	 */
	std::vector<std::unordered_map<node, index>> getEgoNetPartitions();

	std::vector<Graph> getEgoNets();

	void setParameters(std::map<std::string, std::string> const &new_parameters);

	void setGroundTruth(const Cover &gt);

private:

	const Graph &G;
	PartitionFunction localClusterAlgo, globalClusterAlgo;
	std::vector<std::unordered_map<node, index>> egoNetPartitions; // for each node: (global node ID, set ID in ego-net)
	std::vector<count> egoNetPartitionCounts; // number of partitions in the ego-net
	std::vector<node> personaOffsets; // personas of node u are the nodes from [u] to [u+1]-1
	Graph personaGraph; // graph with the split personas
	Partition personaPartition;
	Cover cover; // the result of the algorithm
	mutable std::unordered_map<std::string, double> executionInfo;
	std::vector<Graph> egoNets;
	std::unordered_map<std::string, std::string> parameters;
	AdjacencyArray directedG;
	Cover groundTruth;

	struct Edge {
		node firstNode;
		node secondNode;
		edgeweight weight;

		Edge(node firstNode, node secondNode, edgeweight weight)
				: firstNode(firstNode), secondNode(secondNode), weight(weight) {}
	};

	std::vector<std::vector<EgoSplitting::Edge>> personaEdges; // for each node: edges between its personas

	void init();

	void createEgoNets();

	void splitIntoPersonas();

	void connectPersonas();

	void createPersonaClustering();

	void createCover();

	/**
	 * Connect the personas of a node. Returns a list of edges between the persona indexes.
	 */
	std::vector<Edge> connectEgoPartitionPersonas(
			const Graph &egoGraph, const Partition &egoPartition) const;

};

struct EgoNetData {
	const Graph &G;
	const AdjacencyArray &directedG;
	const Cover &groundTruth;
	node egoNode;
	Graph &egoGraph;
	NodeMapping &egoMapping;
	const std::unordered_map<std::string, std::string> &parameters;
};

} /* namespace NetworKit */

#endif /* EGOSPLITTING_H */
