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

void NodeMapping::reset() {
	for (node v : localToGlobal) {
		globalToLocal[v] = none;
	}
	localToGlobal.clear();
}

void NodeMapping::reset(index end) {
	for (node v = end; v < localToGlobal.size(); ++v)
		globalToLocal[v] = none;
	localToGlobal.resize(end);
}

} /* namespace NetworKit */
