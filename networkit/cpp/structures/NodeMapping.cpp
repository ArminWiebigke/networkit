/*
 * NodeMapping.cpp
 *
 * Created: 2019-03-28
 * Author: Armin Wiebigke
 */

#include "NodeMapping.h"


namespace NetworKit {

NodeMapping::NodeMapping(const NetworKit::Graph &G)
		: globalToLocal(G.upperNodeIdBound(), none) {}

void NodeMapping::addNode(NetworKit::node u) {
	if (!isMapped(u)) {
		globalToLocal[u] = localToGlobal.size();
		localToGlobal.push_back(u);
	}
}

node NodeMapping::local(NetworKit::node globalNode) const {
	return globalToLocal[globalNode];
}

node NodeMapping::global(NetworKit::node localNode) const {
	return localToGlobal[localNode];
}

bool NodeMapping::isMapped(NetworKit::node u) const {
	return local(u) != none;
}

count NodeMapping::nodeCount() const {
	return localToGlobal.size();
}

std::vector<node> NodeMapping::globalNodes() const {
	return localToGlobal;
}

void NodeMapping::reset() {
	for (node v : localToGlobal) {
		globalToLocal[v] = none;
	}
	localToGlobal.clear();
}

} /* namespace NetworKit */