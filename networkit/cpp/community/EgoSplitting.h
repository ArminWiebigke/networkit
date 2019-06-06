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
#include "CommunityDetectionAlgorithm.h"
#include "../structures/Cover.h"
#include "../structures/AdjacencyArray.h"
#include "../structures/NodeMapping.h"

namespace NetworKit {

/**
 * Ego Splitting is a framework to detect overlapping communities.
 * https://dl.acm.org/citation.cfm?id=3098054
 */
class EgoSplitting : public Algorithm {

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
	             std::function<Partition(Graph &)> clusterAlgo);

	/**
	 * Construct an instance of this algorithm for the input graph.
	 *
	 * @param[in]	G   input graph
	 * @param[in]   localClusterAlgo    algorithm to cluster the ego-net
	 * @param[in]   globalClusterAlgo   algorithm to cluster the persona graph
	 */
	EgoSplitting(const Graph &G,
	             std::function<Partition(Graph &)> localClusterAlgo,
	             std::function<Partition(Graph &)> globalClusterAlgo);

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
	 * Get timings for the parts of the algorithm.
	 * @return A map that maps the timer name to its value.
	 */
	std::map<std::string, double> getTimings();

	/**
	 * Get additional information about the execution.
	 * @return A map that maps a name to a value.
	 */
	std::map<std::string, double> getExecutionInfo();

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
	std::function<Partition(Graph &)> localClusterAlgo, globalClusterAlgo;
	std::vector<std::unordered_map<node, index>> egoNetPartitions; // for each node: <global node ID, set ID in ego-net>
	std::vector<count> egoNetPartitionCounts; // number of partitions in the ego-net
	std::vector<node> personaOffsets; // personas of node u are the nodes from [u] to [u+1]-1
	Graph personaGraph; // graph with the split personas
	Partition personaPartition;
	Cover cover; // the result of the algorithm
	std::map<std::string, double> timings;
	std::map<std::string, double> executionInfo;
	std::vector<Graph> egoNets;
	std::map<std::string, std::string> parameters;
//	Graph edgeScoreGraph;
	AdjacencyArray directedEdges;
	Cover groundTruth;

	void init();

	void createEgoNets();

	void splitIntoPersonas();

	void connectPersonas();

	void createPersonaClustering();

	void createCover();

	void extendEgoNet(Graph &egoGraph, node u, NodeMapping &neighbors, count extendNodeCnt);

	Graph weightedEdgesGraph(Graph const &inputGraph);

	void
	consensusWeighting(Graph &egoGraph, node u, const NodeMapping &nodeMapping,
	                   count extendNodeCnt);

	std::vector<std::pair<node, double>> scoreEdgeCount(node u, const NodeMapping &neighbors);

	std::vector<std::pair<node, double>>
    scoreTriangles(node u, const NodeMapping &neighbors,
                   std::vector<std::set<node>> &triangleEdges,
                   Graph const &egoGraph);

	double normalizeScore(node v, double score) const;

	/**
	 * Remove all edges adjacent to nodes with a low degree
	 */
	void removeLowDegreeNodes(Graph &egoGraph, count minDegree, count directNeighborsCnt) const;

	void removeLowTriangleCntNodes(Graph &egoGraph, count directNeighborsCnt) const;


	/**
	 * Search for triangles and execute function for each found triangle
	 */
	void findTriangles(Graph graph, AdjacencyArray directedGraph,
			std::function<void(node, node, node)> triangleFunc) const;

	Partition createGroundTruthPartition(Graph &egoGraph, NodeMapping &mapping, node egoNode) const;
};

} /* namespace NetworKit */


#endif /* EGOSPLITTING_H */
