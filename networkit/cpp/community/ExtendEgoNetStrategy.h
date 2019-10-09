/*
 * ExtendEgoNetStrategy.h
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#ifndef EXTENDE_EGO_NET_STRATEGY_H
#define EXTENDE_EGO_NET_STRATEGY_H

#include <unordered_map>

#include "../base/Algorithm.h"
#include "../Globals.h"
#include "EgoSplitting.h"

namespace NetworKit {

class ExtendEgoNetStrategy : public Algorithm, public Timings {
public:
	ExtendEgoNetStrategy(const EgoNetData &egoNetData, count maxExtendedNodes,
	                     const Graph &egoGraph, node egoNode);

	virtual std::vector<node> getNodes();

protected:
	const Graph &G;
	const AdjacencyArray &directedG;
	const Graph &egoGraph;
	const NodeMapping &egoMapping;
	node egoNode;
	const std::unordered_map<std::string, std::string> &parameters;
	std::vector<node> result;
	count maxExtendedNodes;
};

} /* namespace NetworKit */

#endif //EXTENDE_EGO_NET_STRATEGY_H
