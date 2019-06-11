/*
 * NodeMapping.h
 *
 * Created: 2019-03-28
 * Author: Armin Wiebigke
 */

#ifndef NODEMAPPING_H
#define NODEMAPPING_H

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Map node IDs between a global graph and a local graph
 */
class NodeMapping {
public:
	NodeMapping() = default;

	explicit NodeMapping(Graph const &G);

	// TODO: Allow arbitary node mappings.
	// TODO: Multiple mappings. Use vector<set<node>> for localToGlobal

	/**
	 * Add a node to the mapping. Does nothing if the node is already mapped.
	 * @param u The node to add.
	 * @return true if the node was added, false if node was already mapped
	 */
	bool addNode(node u);

	void addDummy();

	void addMapping(node global, node local);

	/**
	 * Get the local node from a global node.
	 * @param globalNode The global node.
	 * @return The mapped local node.
	 */
	node local(node globalNode) const;

	/**
	 * Get the global node from a local node.
	 * @param localNode The local node.
	 * @return The mapped global node.
	 */
	node global(node localNode) const;

	/**
	 * Check if a global node is mapped.
	 * @param u The global node.
	 * @return True iff the node is mapped.
	 */
	bool isMapped(node u) const;

	/**
	 * Get the number of mapped nodes.
	 * @return The number of mapped nodes.
	 */
	count nodeCount() const;

	/**
	 * Returns the global node IDs for all mapped nodes.
	 * @return A vector of the global nodes.
	 */
	std::vector<node> globalNodes() const;

	void reset();

private:
	std::vector<node> globalToLocal;
	std::vector<node> localToGlobal;
};

} /* namespace NetworKit */

#endif /* EGOSPLITTING_H */
