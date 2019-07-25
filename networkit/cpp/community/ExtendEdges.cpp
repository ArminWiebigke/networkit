/*
 * ExtendEdges.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#include "ExtendEdges.h"

namespace NetworKit {

ExtendEdges::ExtendEdges(const EgoNetData &egoNetData, count maxCandidates)
		: ExtendScore(egoNetData, maxCandidates), nodeScores(egoNetData.nodeScores) {

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

	if (parameters.at("extendOverDirected") == "Yes") {
		for (node v : directNeighbors)
			directedG.forEdgesOf(v, countEdges);
		if (parameters.at("extendDirectedBack") == "Yes") {
			for (node v : candidates) {
				directedG.forEdgesOf(v, [&](node, node w, edgeweight weight) {
					if (isDirectNeighbor(w))
						nodeScores[v] += weight;
				});
			}
		}
	} else {
		for (node v : directNeighbors) {
			G.forEdgesOf(v, countEdges);
		}
	}
	addTime(timer, "3    Count edges");

	std::vector<node> all_candidates(candidates);
	auto newEnd = std::remove_if(candidates.begin(), candidates.end(), [&](node v) {
		return (isDirectNeighbor(v) || v == egoNode);
	});
	candidates.resize(newEnd - candidates.begin());
	addTime(timer, "5    Remove neighbors as candidates");

	// Calculate score for each candidate
	result.reserve(candidates.size());
	for (node v : candidates) {
		double numEdges = nodeScores[v];
		if (numEdges >= 3) {
			double score = normalizeScore(v, numEdges);
			result.emplace_back(v, score);
		}
	}
	addTime(timer, "7    Calculate score");
	std::sort(result.begin(), result.end(), [](NodeScore a, NodeScore b){
		return a.second > b.second;
	});
	if (result.size() > maxExtendedNodes)
		result.resize(maxExtendedNodes);
	addTime(timer, "9    Take best candidates");

	for (node v: all_candidates)
		nodeScores[v] = 0.0;
	addTime(timer, "a    Reset node scores");
}

double ExtendEdges::normalizeScore(node v, double score) const {
	std::string scoreStrategy = parameters.at("scoreStrategy");
	if (scoreStrategy == "constant")
		return 1.0;
	if (scoreStrategy == "random")
		return Aux::Random::real();
	if (scoreStrategy == "score" || scoreStrategy == "none")
		return score * 1.0;
	if (scoreStrategy == "score_normed")
		return score * 1.0 / G.degree(v);
	if (scoreStrategy == "score^1.5_normed")
		return std::pow(score, 1.5) / G.degree(v);
	if (scoreStrategy == "score^2_normed")
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