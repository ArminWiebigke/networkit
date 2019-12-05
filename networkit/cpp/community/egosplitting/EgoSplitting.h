/*
 * egosplitting.h
 *
 * Created: 2018-12-11
 * Author: Armin Wiebigke
 */

#ifndef EGOSPLITTING_H
#define EGOSPLITTING_H

#include <unordered_map>
#include <functional>
#include <ostream>

#include "../../Globals.h"
#include "../../base/Algorithm.h"
#include "../../structures/Cover.h"
#include "../../structures/LowToHighDirectedGraph.h"
#include "../../structures/NodeMapping.h"
#include "../../auxiliary/Timer.h"
#include "../../auxiliary/Timings.h"
#include "../../structures/MemoizationTable.h"
#include "../cleanup/StochasticDistribution.h"
#include "../cleanup/SignificanceCalculator.h"
#include "../../structures/SparseVector.h"
#include "../../auxiliary/ParallelTimings.h"
#include "../ClusteringFunctionFactory.h"

namespace NetworKit {

/**
 * @ingroup community
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
class EgoSplitting : public Algorithm, public ParallelTimings {
public:
	/**
	 * Construct an instance of this algorithm, using the default clustering algorithms.
	 *
	 * @param[in]	G   input graph
	 */
	explicit EgoSplitting(const Graph &G, bool parallelEgoNetEvaluation = true);

	/**
	 * Construct an instance of this algorithm, using the given clustering algorithm for both the
	 * local and the global clustering.
	 *
	 * @param[in]	G   input graph
	 * @param[in]   clusterAlgo    algorithm to cluster the ego-net and the persona graph
	 */
	EgoSplitting(const Graph &G, bool parallelEgoNetEvaluation, ClusteringFunction clusterAlgo);

	/**
	 * Construct an instance of this algorithm, using the given clustering algorithms for the
	 * local and the global clustering.
	 *
	 * @param[in]	G   input graph
	 * @param[in]   localClusterAlgo    algorithm to cluster the ego-net
	 * @param[in]   globalClusterAlgo   algorithm to cluster the persona graph
	 */
	EgoSplitting(const Graph &G, bool parallelEgoNetEvaluation, ClusteringFunction localClusterAlgo,
	             ClusteringFunction globalClusterAlgo);

	void run() override;

	/**
	 * Returns the detected communities.
	 * @return cover containing all detected communities
	 */
	Cover getCover();

	std::string toString() const override;

	/**
	 * Set parameters of the algorithm.
	 * @param new_parameters map that maps parameter name to parameter value
	 */
	void setParameters(std::map<std::string, std::string> const &new_parameters);

	/**
	 * @a Debugging
	 * Get the partitions of the ego-nets. A ego-net partition maps a node to its partition ID.
	 * @return A vector of the ego-net partitions.
	 */
	std::vector<std::unordered_map<node, index>> getEgoNetPartitions();

	/**
	 * @a Debugging
	 * Get the ego-nets. An ego-net is stored as a list of directed edges.
	 * For each node @a u, the list also contains an edge (u,u), which preserves isolated nodes.
	 * @return map that maps each node to its ego-net
	 */
	std::unordered_map<node, std::vector<WeightedEdge>> getEgoNets();


	/**
	 * @a Debugging
	 * Set the ground truth communities.
	 * @param gt cover containing the ground truth communities
	 */
	void setGroundTruth(const Cover &gt);

private:
	const Graph &G;
	bool parallelEgoNetEvaluation;
	ClusteringFunction localClusteringAlgo, globalClusteringAlgo;
	std::vector<std::unordered_map<node, index>> egoNetPartitions; // for each node: (global node ID, set ID in ego-net)
	// for each node: (global node ID, set ID in ego-net).
	// Includes nodes of the extended ego-net. These partitions are not by the algorithm itself,
	// as we only need the nodes of the original ego-net. The partitions are only useful for
	// the analysis of the algorithm (use getEgoNetPartitions() to get them).
	std::vector<std::unordered_map<node, index>> egoNetExtendedPartitions;
	std::vector<count> egoNetPartitionCounts; // number of partitions in the ego-net
	std::vector<node> personaOffsets; // personas of node u are the nodes from [u] to [u+1]-1
	Graph personaGraph;
	Partition personaPartition;
	Cover resultCover; // the result of the algorithm
	std::unordered_map<node, std::vector<WeightedEdge>> egoNets;
	std::unordered_map<std::string, std::string> parameters;
	LowToHighDirectedGraph directedG;
	Cover groundTruth;
	std::vector<std::vector<WeightedEdge>> personaEdges; // for each node: edges between its personas
	StochasticDistribution stochasticDistribution; // Stochastic distribution needed for significance calculations

	void init();

	void createEgoNets();
	void createEgoNetsSimple();

	void splitIntoPersonas();

	void connectPersonas();

	void createPersonaClustering();

	std::vector<std::vector<node>> getCommunitiesFromPersonaClustering();

	void createCover(const std::vector<std::vector<node>> &communities);

	/**
	 * Connect the personas of a node. Returns a list of edges between the persona indexes.
	 */
	std::vector<WeightedEdge> connectEgoPartitionPersonas(
			const Graph &egoGraph, const Partition &egoPartition) const;

	void storeEgoNet(const Graph &egoGraph, const NodeMapping &egoMapping, node egoNode);

	void cleanUpCommunities(std::vector<std::vector<node>> &communities);

	void discardSmallCommunities(std::vector<std::vector<node>> &communities);
};

/**
 * Contains data that is used for the execution of the algorithm.
 * For performance reasons, most of the data is persistent across the entire algorithm.
 */
struct EgoNetData {
	// Global data
	const Graph &G;
	const LowToHighDirectedGraph &directedG;
	const Cover &groundTruth;
	const std::unordered_map<std::string, std::string> &parameters;
	MemoizationTable<double> &sigTable;
	// Ego-net local data
	NodeMapping &egoMapping;
	SparseVector<double> &nodeScores;
	SparseVector<node> &significantGroup;
	SparseVector<std::vector<count>> &edgesToGroups;
	SignificanceCalculator &significanceCalculator;
};

} /* namespace NetworKit */

#endif /* EGOSPLITTING_H */
