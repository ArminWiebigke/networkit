/*
 * ExtendScore.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#include "ExtendScore.h"

namespace NetworKit {

ExtendScore::ExtendScore(const EgoNetData &egoNetData)
		: G(egoNetData.G),
		  directedG(egoNetData.directedG),
		  egoGraph(egoNetData.egoGraph),
		  egoMapping(egoNetData.egoMapping),
		  egoNode(egoNetData.egoNode),
		  parameters(egoNetData.parameters) {

}

std::vector<std::pair<node, double>> ExtendScore::getScores() {
	hasRun = false;
	return std::move(result);
}

} /* namespace NetworKit */
