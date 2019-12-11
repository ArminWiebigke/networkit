/*
 * EgoSplitting.h
 *
 * Created: 2019-05-28
 * Author: Armin Wiebigke
 */

#ifndef SLPA_H
#define SLPA_H

#include <unordered_map>
#include <string>

#include "../graph/Graph.h"
#include "../Globals.h"
#include "CommunityDetectionAlgorithm.h"
#include "../structures/Partition.h"
#include "../structures/Cover.h"

namespace NetworKit {
using label = index;

using LabelCounts = std::unordered_map<label, count>;

/**
 * A Speaker-listener Label Propagation Algorithm (SLPA) for overlapping community detection.
 * Each node has a memory and remembers all labels it has ever seen. In one iteration, all nodes
 * are selected in a random order. Each node receives one label from each of its neighbors.
 * The label that a node sends is randomly selected from the labels in its memory. The chance to
 * select a label is proportional to the number of times the label is in the memory. Then one
 * of the received labels is added to the memory of the node.
 * In the end, each node is assigned to the communities of the labels that exceed a given threshold.
 // TODO: cite
 */
class SLPA : public Algorithm {

public:
	/**
	 * Construct an instance of this algorithm for the input graph.
	 * @param graph Input graph
	 * @param threshold Threshold to accept labels
	 * @param numIterations Number of iterations
	 */
	explicit SLPA(const Graph &graph, double threshold = 0.1, count numIterations = 100);

	/**
	 * Construct an instance of this algorithm for the input graph, based on an existing partition.
	 * @param graph Input graph
	 * @param basePartition Base partition
	 * @param threshold Threshold to accept labels
	 * @param numIterations Number of iterations
	 */
	SLPA(const Graph &graph, const Partition &basePartition, double threshold = 0.1,
	     count numIterations = 100);

	void run() override;

	std::string toString() const override;

	bool isParallel() const override;


	/**
	 * Returns a cover of the graph. Each node is assigned all labels that exceed the threshold.
	 * @return The cover
	 */
	Cover getCover();

	/**
	 * Returns a partition of the graph. Each nodes is assigned the most frequently heard label.
	 * @return The partition
	 */
	Partition getPartition();

private:
	const Graph &graph; // Input graph
	count numIterations; // Number of iteratitions
	double threshold; // Threshold to accept a label
	Cover cover; // Result cover
	Partition partition; // Result partition

	/**
	 * For one node, get a label from each neighbor.
	 * @param u The node that is listening.
	 * @return The received labels and their number of occurrences.
	 */
	LabelCounts listen(node u);

	/**
	 * Select the label with the most occurrences.
	 * @param labelCounts The labels and their number of occurrences.
	 * @return The best label
	 */
	label selectLabel(const LabelCounts& labelCounts);

	/**
	 * Create the result cover and partition.
	 */
	void createResult();

	/**
	 * Send a label from a node.
	 * @param u The node that sends a label.
	 * @return A label
	 */
	label sendLabel(node u) const;

	/**
	 * Add a label to the memory of a node.
	 * @param u The node that receives the label.
	 * @param l The label to add.
	 */
	void addLabelTo(node u, label l);

	std::vector<LabelCounts> nodesMemory;
};


} /* namespace NetworKit */


#endif /* SLPA_H */
