/*
 * ExtendEgoNetStrategy.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#include "ExtendEgoNetStrategy.h"

namespace NetworKit {

ExtendEgoNetStrategy::ExtendEgoNetStrategy(const EgoNetData &egoNetData, count maxExtendedNodes,
                                           const Graph &egoGraph, node egoNode)
		: G(egoNetData.G),
		  directedG(egoNetData.directedG),
		  egoGraph(egoGraph),
		  egoMapping(egoNetData.egoMapping),
		  egoNode(egoNode),
		  parameters(egoNetData.parameters),
		  maxExtendedNodes(maxExtendedNodes) {

}

std::vector<node> ExtendEgoNetStrategy::getNodes() {
	if (!hasRun)
		throw std::runtime_error("Run the algorithm first!");
	hasRun = false;
	return std::move(significantCandidates);
}

} /* namespace NetworKit */
