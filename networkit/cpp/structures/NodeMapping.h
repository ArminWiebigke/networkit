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

	node local(node globalNode);

	node global(node localNode);

	bool isMapped(node u);

	count nodeCount();

	std::vector<node> globalNodes();

	void reset();

private:
	std::vector<node> globalToLocal;
	std::vector<node> localToGlobal;
};

} /* namespace NetworKit */

#endif /* EGOSPLITTING_H */
