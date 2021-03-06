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
	bool addNode(node u) {
		if (!isMapped(u)) {
			globalToLocal[u] = localToGlobal.size();
			localToGlobal.push_back(u);
			return true;
		}
		return false;
	}

	void addDummy() {
		localToGlobal.push_back(none);
	}

	void addMapping(node global, node local) {
		globalToLocal[global] = local;
		localToGlobal[local] = global;
	}

	/**
	 * Get the local node from a global node.
	 * @param globalNode The global node.
	 * @return The mapped local node.
	 */
	node toLocal(node globalNode) const {
		assert(globalToLocal[globalNode] != none);
		return globalToLocal[globalNode];
        }

	/**
	 * Get the global node from a local node.
	 * @param localNode The local node.
	 * @return The mapped global node.
	 */
	node toGlobal(node localNode) const {
		assert(localToGlobal[localNode] != none);
		return localToGlobal[localNode];
	}

	/**
	 * Check if a global node is mapped.
	 * @param globalNode The global node.
	 * @return True iff the node is mapped.
	 */
	bool isMapped(node globalNode) const {
		return globalToLocal[globalNode] != none;
	}

	/**
	 * Get the number of mapped nodes.
	 * @return The number of mapped nodes.
	 */
	count nodeCount() const {
		return localToGlobal.size();
	}

	/**
	 * Returns the global node IDs for all mapped nodes.
	 * @return A vector of the global nodes.
	 */
	const std::vector<node>& globalNodes() const {
		return localToGlobal;
	}

	void reset();

	/**
	 * Reset the mapping while keeping only the smallest local nodes.
	 * @param end The first local node that will be deleted.
	 */
	void reset(index end);

private:
	std::vector<node> globalToLocal;
	std::vector<node> localToGlobal;
};

} /* namespace NetworKit */

#endif /* EGOSPLITTING_H */
