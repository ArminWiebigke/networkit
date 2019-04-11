//
// Created by armin on 3/28/19.
//

#ifndef NODEMAPPING_H
#define NODEMAPPING_H

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Map node IDs between a global graph and a local graph
 */
class NodeMapping {
public:
	explicit NodeMapping(Graph const &G);

	void addNode(node u);

	node local(node globalNode) const;

	node global(node localNode) const;

	bool isMapped(node u) const;

	count nodeCount() const;

	std::vector<node> globalNodes() const;

	void reset();

private:
	std::vector<node> globalToLocal;
	std::vector<node> localToGlobal;
};

} /* namespace NetworKit */

#endif /* EGOSPLITTING_H */
