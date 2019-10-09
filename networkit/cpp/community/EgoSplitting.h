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
#include "../structures/MemoizationTable.h"

namespace NetworKit {

/**
 * Ego Splitting is a framework to detect overlapping communities.
 * The ego-net of each node is partitioned by a clustering algorithm. The ego-net of a node u
 * (u is the ego-node) is the subgraph induced by the neighbors of u. For each detected subset of
 * the ego-net, a copy of the ego-node is created, a so-called persona.
 * After analyzing the ego-nets of all nodes, a persona graph is created that consists of the
 * personas of all nodes. Each edge in the input graph corresponds to exactly one edge in the
 * persona graph. A second clustering algorithm is used on the persona graph to detect
 * non-overlapping communities. Each node is then assigned all communities of its personas,
 * resulting in (possibly) overlapping communites.
 * https://dl.acm.org/citation.cfm?id=3098054
 */
class EgoSplitting : public Algorithm, public Timings {
	using PartitionFunction = std::function<Partition(const Graph &)>;

public:
	/**
	 * Construct an instance of this algorithm, using the default clustering algorithms.
	 *
	 * @param[in]	G   input graph
	 */
	explicit EgoSplitting(const Graph &G);

	/**
	 * Construct an instance of this algorithm, using the given clustering algorithm for both the
	 * local and the global clustering.
	 *
	 * @param[in]	G   input graph
	 * @param[in]   clusteringAlgo    algorithm to cluster the ego-net and the persona graph
	 */
	EgoSplitting(const Graph &G,
	             PartitionFunction clusteringAlgo);

	/**
	 * Construct an instance of this algorithm, using the given clustering algorithms for the
	 * local and the global clustering.
	 *
	 * @param[in]	G   input graph
	 * @param[in]   localClusteringAlgo    algorithm to cluster the ego-net
	 * @param[in]   globalClusteringAlgo   algorithm to cluster the persona graph
	 */
	EgoSplitting(const Graph &G,
	             PartitionFunction localClusteringAlgo,
	             PartitionFunction globalClusteringAlgo);

	void run() override;

	/**
	 * Returns the detected communities.
	 * @return cover containing all detected communities
	 */
	Cover getCover();

	std::string toString() const override;

	/**
	 * Get the partitions of the ego-nets. A ego-net partition maps a node to its partition ID.
	 * @return A vector of the ego-net partitions.
	 */
	std::vector<std::unordered_map<node, index>> getEgoNetPartitions();

	std::unordered_map<node, std::vector<WeightedEdge>> getEgoNets();

	void setParameters(std::map<std::string, std::string> const &new_parameters);

	void setGroundTruth(const Cover &gt);

private:

	const Graph &G;
	PartitionFunction localClusteringAlgo, globalClusteringAlgo;
	std::vector<std::unordered_map<node, index>> egoNetPartitions; // for each node: (global node ID, set ID in ego-net)
	// for each node: (global node ID, set ID in ego-net).
	// Includes nodes of the extended ego-net. These partitions are not by the algorithm itself,
	// as we only need the nodes of the original ego-net. The partitions are only useful for
	// the analysis of the algorithm (use getEgoNetPartitions() to get them).
	std::vector<std::unordered_map<node, index>> egoNetExtendedPartitions;
	std::vector<count> egoNetPartitionCounts; // number of partitions in the ego-net
	std::vector<node> personaOffsets; // personas of node u are the nodes from [u] to [u+1]-1
	Graph personaGraph; // graph with the split personas
	Partition personaPartition;
	Cover cover; // the result of the algorithm
	std::unordered_map<node, std::vector<WeightedEdge>> egoNets;
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

	void storeEgoNet(const Graph &egoGraph, const NodeMapping &egoMapping, node egoNode);

};

/**
 * Contains data that is used for the execution of the algorithm.
 * For performance reasons, most of the data is persistent across the entire algorithm.
 */
struct EgoNetData {
	const Graph &G;
	const AdjacencyArray &directedG;
	const Cover &groundTruth;
	NodeMapping &egoMapping;
	const std::unordered_map<std::string, std::string> &parameters;
	MemoizationTable<double> &sigTable;
	std::vector<double> &nodeScores;
	std::vector<node> &significantGroup;
	std::vector<std::vector<count>> &edgesToGroups;
};

} /* namespace NetworKit */

#endif /* EGOSPLITTING_H */
