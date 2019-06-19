/*
 * ExtendEdges.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#include "ExtendEdges.h"

namespace NetworKit {

ExtendEdges::ExtendEdges(const EgoNetData &egoNetData)
		: ExtendScore(egoNetData) {

}

void ExtendEdges::run() {
	Aux::Timer timer;
	timer.start();
	std::vector<double> edgeScores(G.upperNodeIdBound()); //TODO: Don't allocate this every time
	addTime(timer, "331    Create vector");
	std::vector<node> candidates;
	// Search for all edges to neighbors of neighbors
	std::vector<node> directNeighbors = egoMapping.globalNodes();
	auto isDirectNeighbor = [&](node x) {
		return egoMapping.isMapped(x);
	};
	auto countEdges = [&](node v, node w, edgeweight weight) {
		if (edgeScores[w] == 0.0)
			candidates.push_back(w);
		edgeScores[w] += weight;
	};

	if (parameters.at("extendOverDirected") == "Yes") {
		for (node v : directNeighbors)
			directedG.forEdgesOf(v, countEdges);
		if (parameters.at("extendDirectedBack") == "Yes") {
			for (node v : candidates) {
				directedG.forEdgesOf(v, [&](node, node w, edgeweight weight) {
					if (isDirectNeighbor(w))
						edgeScores[v] += weight;
				});
			}
		}
	} else {
		for (node v : directNeighbors) {
			G.forEdgesOf(v, countEdges);
		}
	}
	addTime(timer, "333    Count edges");

	auto newEnd = std::remove_if(candidates.begin(), candidates.end(), [&](node v) {
		return (isDirectNeighbor(v) || v == egoNode);
	});
	candidates.resize(newEnd - candidates.begin());
	addTime(timer, "335    Remove neighbors as candidates");

	// Calculate score for each candidate
	result.reserve(candidates.size());
	for (node v : candidates) {
		double numEdges = edgeScores[v];
		double score = normalizeScore(v, numEdges);
		if (numEdges < 3)
			score = 0.0;
		result.emplace_back(v, score);
	}
	addTime(timer, "337    Calculate score");
}

double ExtendEdges::normalizeScore(node v, double score) const {
	std::string scoreStrategy = parameters.at("scoreStrategy");
	if (scoreStrategy == "constant")
		return 1.0;
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

} /* namespace NetworKit */