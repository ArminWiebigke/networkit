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

#include "../Globals.h"
#include "CommunityDetectionAlgorithm.h"
#include "../structures/Cover.h"

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
	 * Construct an instance of this algorithm for the input graph.
	 *
	 * @param[in]	G   input graph
	 * @param[in]   localClusterAlgo    algorithm to cluster the ego-net
	 * @param[in]   globalClusterAlgo   algorithm to cluster the persona graph
	 */
	EgoSplitting(const Graph &G,
				 std::function<Partition(Graph &)> localClusterAlgo,
				 std::function<Partition(Graph &)> globalClusterAlgo,
				 const Cover &groundTruth);

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

	std::map<std::string, double> getTimings();

	std::map<std::string, double> getExecutionInfo();

	std::vector<std::unordered_map<node, index>> getEgoNetPartitions();

private:

	const Graph &G;
	std::function<Partition(Graph &)> localClusterAlgo, globalClusterAlgo;
	std::vector<std::unordered_map<node, index>> egoNetPartitions; // for each node: <global node ID, set ID in ego-net>
	std::vector<node> personaOffsets; // personas of node u are the nodes from [u] to [u+1]-1
	Graph personaGraph; // graph with the split personas
	Partition personaPartition;
	Cover cover; // the result of the algorithm
	std::map<std::string, double> timings;
	std::map<std::string, double> executionInfo;
	Cover groundTruth;
	std::vector<Partition> egoPartitions;

	void init();

	void createEgoNets();

	void splitIntoPersonas();

	void connectPersonas();

	void createPersonaClustering();

	void createCover();
};

} /* namespace NetworKit */


#endif /* EGOSPLITTING_H */