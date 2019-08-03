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

	if (parameters.at("Edges Iterative") == "Yes") {
		count iterations = 10;
		count maxNodesAdded = maxExtendedNodes / iterations + 1;
		std::vector<std::pair<double, node>> candidateScores;
		for (node v : candidates) {
			double numEdges = nodeScores[v];
			if (numEdges >= 3) {
				double score = normalizeScore(v, numEdges);
				candidateScores.emplace_back(score, v);
			}
		}

		for (count i = 0; i < iterations; ++i) {
			std::sort(candidateScores.begin(), candidateScores.end());
			std::vector<node> addedCandidates;
			for (count j = 0; j < maxNodesAdded; ++j) {
				if (candidateScores.empty())
					break;
				auto pair = candidateScores.back();
				node v = pair.second;
				result.emplace_back(v, pair.first);
				nodeScores[v] = 0;
				addedCandidates.push_back(v);
				candidateScores.pop_back();
			}
			if (candidateScores.empty())
				break;

			// Update number of edges
			for (node v : addedCandidates) {
				G.forEdgesOf(v, [&](node, node w, edgeweight weight) {
					if (nodeScores[w] > 0)
						nodeScores[w] += weight;
				});
			}

			// Update remaining candidates
			for (auto &pair : candidateScores) {
				node v = pair.second;
				double numEdges = nodeScores[v];
				pair.first = normalizeScore(v, numEdges);
			}
		}
	} else {
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
	}

	std::sort(result.begin(), result.end(), [](NodeScore a, NodeScore b) {
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
	std::string scoreStrategy = parameters.at("Edges Score Strategy");
	if (scoreStrategy == "constant")
		return 1.0;
	if (scoreStrategy == "Random")
		return Aux::Random::real();
	if (scoreStrategy == "Edges" || scoreStrategy == "none")
		return score * 1.0;
	if (scoreStrategy == "Edges / Degree")
		return score * 1.0 / G.degree(v);
	if (scoreStrategy == "Edges^2 / Degree")
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