/*
 * ExtendEdges.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#include "ExtendEdges.h"

namespace NetworKit {

ExtendEdges::ExtendEdges(const EgoNetData &egoNetData, count maxCandidates,
                         const Graph &egoGraph, node egoNode)
		: ExtendEgoNetStrategy(egoNetData, maxCandidates, egoGraph, egoNode),
		  nodeScores(egoNetData.nodeScores) {

}

void ExtendEdges::run() {
	Aux::Timer timer;
	timer.start();
	std::vector<node> candidates;
	// Search for all edges to neighbors of neighbors
	std::vector<node> directNeighbors = egoMapping.globalNodes();
	auto isDirectNeighbor = [&](node x) {
		return egoMapping.isMapped(x);
	};
	auto countEdges = [&](node v, node w, edgeweight weight) {
		if (nodeScores[w] == 0.0)
			candidates.push_back(w);
		nodeScores[w] += weight;
	};
	for (node v : directNeighbors) {
		G.forEdgesOf(v, countEdges);
	}
	addTime(timer, "3    Count edges");

	std::vector<node> all_candidates(candidates);
	auto shouldBeRemoved = [&](node v){
		return isDirectNeighbor(v) || v == egoNode;
	};
	auto newEnd = std::remove_if(candidates.begin(), candidates.end(), [&](node v) {
		return shouldBeRemoved(v);
	});
	candidates.resize(newEnd - candidates.begin());
	addTime(timer, "5    Remove neighbors as candidates");


	// Calculate score for each candidate
	using NodeAndScore = std::pair<node, double>;
	std::vector<NodeAndScore> candidatesAndScores;
	candidatesAndScores.reserve(candidates.size());
	for (node v : candidates) {
		double numEdges = nodeScores[v];
		if (numEdges >= 3) {
			double score = normalizeScore(v, numEdges);
			candidatesAndScores.emplace_back(v, score);
		}
	}
	addTime(timer, "7    Calculate score");


	std::sort(candidatesAndScores.begin(), candidatesAndScores.end(),
	          [](NodeAndScore a, NodeAndScore b) {
		          return a.second > b.second;
	          });
	if (candidatesAndScores.size() > maxExtendedNodes)
		candidatesAndScores.resize(maxExtendedNodes);
	for (NodeAndScore nodeAndScore : candidatesAndScores) {
		result.push_back(nodeAndScore.first);
	}
	addTime(timer, "9    Take best candidates");

	for (node v: all_candidates)
		nodeScores[v] = 0.0;
	addTime(timer, "a    Reset node scores");

	hasRun = true;
}

double ExtendEdges::normalizeScore(node v, double score) const {
	std::string scoreStrategy = parameters.at("Edges Score Strategy");
	if (scoreStrategy == "constant")
		return 1.0;
	if (scoreStrategy == "Random")
		return Aux::Random::real();
	if (scoreStrategy == "Edges" || scoreStrategy == "none")
		return score * 1.0;
	if (scoreStrategy == "Edges div Degree")
		return score * 1.0 / G.degree(v);
	if (scoreStrategy == "Edges pow 2 div Degree")
		return score * score * 1.0 / G.degree(v);
	throw std::runtime_error(scoreStrategy + " is not a valid score strategy!");
}

bool ExtendEdges::isParallel() const {
	return false;
}

std::string ExtendEdges::toString() const {
	return "ExtendEdges";
}

} /* namespace NetworKit */