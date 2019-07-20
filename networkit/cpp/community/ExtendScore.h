/*
 * ExtendScore.h
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#ifndef EXTENDSCORE_H
#define EXTENDSCORE_H

#include <unordered_map>

#include "../base/Algorithm.h"
#include "../Globals.h"
#include "EgoSplitting.h"

namespace NetworKit {

class ExtendScore : public Algorithm, public Timings {
public:
	using NodeScore = std::pair<node, double>;
	ExtendScore(const EgoNetData &egoNetData, count maxExtendedNodes);

	virtual std::vector<NodeScore> getScores();

protected:
	const Graph &G;
	const AdjacencyArray &directedG;
	const Graph &egoGraph;
	const NodeMapping &egoMapping;
	node egoNode;
	const std::unordered_map<std::string, std::string> &parameters;
	std::vector<NodeScore> result;
	count maxExtendedNodes;
};

} /* namespace NetworKit */

#endif //EXTENDSCORE_H
