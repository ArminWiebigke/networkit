/*
 * LowToHighDirectedGraph.h
 *
 * Created: 2018-12-12
 * Author: Armin Wiebigke
 */

#ifndef LOW_TO_HIGH_DIRECTED_GRAPH_H
#define LOW_TO_HIGH_DIRECTED_GRAPH_H

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * This class takes an undirected graph and transforms it into a directed graph.
 * Each undirected edge {u,v} is transformed into a single directed edge (u,v) from the lower degree
 * node to the higher degree node, i.e. degree(u) <= degree(v).
 * On the directed graph, some algorithms, e.g. triangle counting, can be performed much faster
 * than on the undirected graph.
 */
class LowToHighDirectedGraph {
public:
	explicit LowToHighDirectedGraph(const Graph &G);

	LowToHighDirectedGraph() = default;

	/**
	 * Iterate over all incident edges of a node and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameters <code>(node, node, edgeweight)</code>,  where the first node
	 * is @a u and the second is a neighbor of @a u.
	 */
	template<typename L>
	void forEdgesOf(node u, L handle) const;

private:
	// The directed graph is stored as an adjacency array.
	std::vector<index> edgesBegin;
	std::vector<node> edgeTargets;
	std::vector<edgeweight> edgeWeights;
};

template<typename L>
void LowToHighDirectedGraph::forEdgesOf(node u, L handle) const {
	for (index i = edgesBegin[u]; i < edgesBegin[u + 1]; ++i) {
		handle(u, edgeTargets[i], edgeWeights[i]);
	}
}

} /* namespace NetworKit */

#endif //LOW_TO_HIGH_DIRECTED_GRAPH_H
