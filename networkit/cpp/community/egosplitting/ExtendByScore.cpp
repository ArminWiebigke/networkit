/*
 * ExtendByScore.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#include "ExtendByScore.h"

namespace NetworKit {

ExtendByScore::ExtendByScore(const EgoNetData &egoNetData, count maxCandidates,
                             const Graph &egoGraph, node egoNode)
		: ExtendEgoNetStrategy(egoNetData, maxCandidates, egoGraph, egoNode),
		  nodeScores(egoNetData.nodeScores), scoreStrategy(parameters.at("Edges Score Strategy")),
		  significance(egoNetData.stochasticSignificance) {

}

void ExtendByScore::run() {
	Aux::Timer timer;
	timer.start();

	std::vector<node> candidates = searchForCandidates();
	addTime(timer, "3    Count edges");

	std::vector<NodeAndScore> candidatesAndScores = calculateScores(candidates);
	addTime(timer, "7    Calculate score");

	takeBestCandidates(candidatesAndScores);
	addTime(timer, "9    Take best candidates");

	nodeScores.reset();
	addTime(timer, "a    Reset node scores");

	hasRun = true;
}

std::vector<node> ExtendByScore::searchForCandidates() {
	std::vector<node> candidates;
	std::vector<node> egoNetNodes = egoMapping.globalNodes();
	auto isInEgoNet = [&](node x) {
		return egoMapping.isMapped(x);
	};
	count internalStubs = 0;
	outgoingStubs = 0;
	auto countEdges = [&](node egoNetNode, node neighbor, edgeweight weight) {
		if (!nodeScores.indexIsUsed(neighbor)) {
			candidates.push_back(neighbor);
			nodeScores.insert(neighbor, 0.0);
		}
		nodeScores[neighbor] += weight;
		if (isInEgoNet(neighbor)) {
			++internalStubs;
		} else {
			++outgoingStubs;
		}
	};
	for (node egoNetNode : egoNetNodes) {
		G.forEdgesOf(egoNetNode, countEdges);
	}
	externalStubs = G.numberOfEdges() * 2 - internalStubs - outgoingStubs;

	// Remove ego-net nodes and ego-node as candidates
	auto shouldBeRemoved = [&](node v) {
		return isInEgoNet(v) || v == egoNode;
	};
	auto newEnd = std::remove_if(candidates.begin(), candidates.end(), [&](node v) {
		return shouldBeRemoved(v);
	});
	candidates.resize(newEnd - candidates.begin());

	return candidates;
}

std::vector<ExtendByScore::NodeAndScore>
ExtendByScore::calculateScores(const std::vector<node> &candidates) const {
	std::vector<NodeAndScore> candidatesAndScores;
	candidatesAndScores.reserve(candidates.size());
	for (node candidate : candidates) {
		double numEdges = nodeScores[candidate];
		if (numEdges >= 3) {
			double score = calculateScore(candidate, numEdges);
			candidatesAndScores.emplace_back(candidate, score);
		}
	}
	return candidatesAndScores;
}

double ExtendByScore::calculateScore(node v, count numEdges) const {
	if (scoreStrategy == "constant")
		return 1.0;
	if (scoreStrategy == "Random")
		return Aux::Random::real();
	if (scoreStrategy == "Edges" || scoreStrategy == "none")
		return (double) numEdges;
	if (scoreStrategy == "Edges div Degree")
		return (double) numEdges / G.degree(v);
	if (scoreStrategy == "Edges pow 2 div Degree")
		return (double) numEdges * numEdges / G.degree(v);
	if (scoreStrategy == "Significance") {
		double rScore = significance.rScore(G.degree(v), numEdges, outgoingStubs, externalStubs);
		return -rScore; // low r-score is better
	}
	throw std::runtime_error(scoreStrategy + " is not a valid score strategy!");
}

void ExtendByScore::takeBestCandidates(std::vector<NodeAndScore> &candidatesAndScores) {
	std::sort(candidatesAndScores.begin(), candidatesAndScores.end(),
	          [](NodeAndScore a, NodeAndScore b) {
		          return a.second > b.second;
	          });
	if (candidatesAndScores.size() > maxExtendedNodes)
		candidatesAndScores.resize(maxExtendedNodes);
	for (NodeAndScore nodeAndScore : candidatesAndScores) {
		significantCandidates.push_back(nodeAndScore.first);
	}
}

bool ExtendByScore::isParallel() const {
	return false;
}

std::string ExtendByScore::toString() const {
	return "ExtendByScore";
}

} /* namespace NetworKit */