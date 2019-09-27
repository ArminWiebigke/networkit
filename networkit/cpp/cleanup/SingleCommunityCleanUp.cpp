/*
 * SingleCommunityCleanup.cpp
 *
 * Created: 2019-09-26
 * Author: Armin Wiebigke
 */
#include "SingleCommunityCleanUp.h"

namespace NetworKit {

using Community = SingleCommunityCleanUp::Community;

Community
SingleCommunityCleanUp::clean(const Community &inputCommunity) {
	Community firstPhaseResult = calculateSignificantNodes(inputCommunity, true);
	Community cleanedCommunity = calculateSignificantNodes(firstPhaseResult, false);
	bool changedDrastically = smallOverlap(inputCommunity, cleanedCommunity);
	if (changedDrastically)
		cleanedCommunity = {};
	return cleanedCommunity;
}

bool SingleCommunityCleanUp::smallOverlap(const Community &inputCommunity,
                                          const Community &cleanedCommunity) const {
	Community intersection;
	std::set_intersection(inputCommunity.begin(), inputCommunity.end(),
	                      cleanedCommunity.begin(), cleanedCommunity.end(),
	                      std::inserter(intersection, intersection.begin()));
	count largerSize = std::max(inputCommunity.size(), cleanedCommunity.size());
	double overlapRatio = (double) intersection.size() / largerSize;
//	std::cout << inputCommunity.size() << " -> " << cleanedCommunity.size()
//	          << (overlapRatio < minOverlapRatio ? " -> bad" : "") << std::endl;
	return overlapRatio < minOverlapRatio;
}

SingleCommunityCleanUp::SingleCommunityCleanUp(const Graph &graph, double scoreThreshold,
                                               double significanceThreshold, double minOverlapRatio)
		: graph(graph),
		  scoreThreshold(scoreThreshold),
		  significanceThreshold(significanceThreshold),
		  minOverlapRatio(minOverlapRatio),
		  edgesToCommunity(graph.upperNodeIdBound()),
		  isInCommunity(graph.upperNodeIdBound()),
		  isCandidate(graph.upperNodeIdBound()),
		  stochastic(2 * graph.numberOfEdges() + graph.numberOfNodes()) {
}

Community
SingleCommunityCleanUp::calculateSignificantNodes(
		const Community &inputCommunity, bool includeNeighbors) {
	community = inputCommunity;
	auto candidates = getCandidatesAndSetUpCalculation(includeNeighbors);
	Community cleanedCommunity;
	while (!community.empty()) {
		auto candidateScores = calculateCandidateScores(candidates);
		cleanedCommunity = findSignificantCandidates(candidateScores);
		bool communityIsSignificant = !cleanedCommunity.empty();
		if (communityIsSignificant)
			break;
		removeWorstNode(candidateScores);
	}
	reset(candidates);
	return cleanedCommunity;
}

std::vector<node>
SingleCommunityCleanUp::getCandidatesAndSetUpCalculation(bool includeNeighbors) {
	std::vector<node> candidates;
	auto addCandidate = [&](node u) {
		candidates.push_back(u);
		isCandidate[u] = true;
	};
	for (node u : community) {
		isInCommunity[u] = true;
		addCandidate(u);
	}
	outgoingCommunityStubs = 0;
	totalCommunityStubs = 0;
	for (node u : community) {
		graph.forNeighborsOf(u, [&](node neighbor) {
			totalCommunityStubs += 1;
			if (!isInCommunity[neighbor])
				outgoingCommunityStubs += 1;
			if (isInCommunity[neighbor] || includeNeighbors) {
				if (!isCandidate[neighbor]) {
					assert(edgesToCommunity[neighbor] == 0);
					addCandidate(neighbor);
				}
				edgesToCommunity[neighbor] += 1;
			}
		});
	}
	// TODO: Is the community part of external nodes? else we have more candidates than external nodes, which causes the order statistics to fail
//	externalNodes = graph.numberOfNodes() - community.size();
	externalNodes = graph.numberOfNodes();

	externalStubs = 2 * graph.numberOfEdges() - totalCommunityStubs;
	return candidates;
}

void
SingleCommunityCleanUp::reset(const std::vector<node> &candidates) {
	for (node u : candidates) {
		edgesToCommunity[u] = 0;
		isCandidate[u] = 0;
	}
	for (node u : community) {
		isInCommunity[u] = false;
	}
	community.clear();
	outgoingCommunityStubs = 0;
	totalCommunityStubs = 0;
	externalNodes = 0;
	externalStubs = 0;
}

// remove the node with the worst (= highest) score from the community
void SingleCommunityCleanUp::removeWorstNode(
		std::vector<ScoreStruct> scores) {
	auto worstNodeIt = std::find_if(scores.rbegin(), scores.rend(),
	                                [&](const ScoreStruct &s) { return isInCommunity[s.candidate]; });
	node worstNode = worstNodeIt->candidate;
	assert(isInCommunity[worstNode]);
	community.erase(worstNode);
	isInCommunity[worstNode] = false;
//	 externalNodes += 1; // TODO: Do this or not?
	count degree = graph.degree(worstNode);
	outgoingCommunityStubs += 2 * edgesToCommunity[worstNode] - degree;
	totalCommunityStubs -= degree;
	externalStubs += degree;
	graph.forNeighborsOf(worstNode, [&](node u) {
		if (isCandidate[u]) {
			assert(edgesToCommunity[u] > 0);
			edgesToCommunity[u] -= 1;
		}
	});
}

Community
SingleCommunityCleanUp::findSignificantCandidates(std::vector<ScoreStruct> scores) const {
	int position = 1;
	int significantNodesCount = 0;
	for (auto scoreStruct : scores) {
		double score = scoreStruct.sScore;
		// significance is the probability Omega_{position}(score, externalNodes)
		double significance = stochastic.orderStatistic(score, externalNodes, position);
		if (significance < significanceThreshold) {
			// The score is much better than expected in the null model, so the node is significant
			significantNodesCount = position;
		} else {
			if (significantNodesCount != 0) {
				// Stop after finding significant nodes and then not significant ones
				break;
			}
		}
		++position;
	}

	Community significantNodes;
	for (count i = 0; i < significantNodesCount; ++i) {
		significantNodes.insert(scores[i].candidate);
	}
	return significantNodes;
}

std::vector<SingleCommunityCleanUp::ScoreStruct>
SingleCommunityCleanUp::calculateCandidateScores(const std::vector<node> &candidates) const {
	std::vector<ScoreStruct> scores;
	for (node u : candidates) {
		// Calculate scores as if the node is not part of the community
		count degree = graph.degree(u);
		count adjustedExternalStubs = externalStubs;
		count adjustedCommOut = outgoingCommunityStubs;
		if (isInCommunity[u]) {
			adjustedExternalStubs += degree;
			adjustedCommOut += 2 * edgesToCommunity[u] - degree;
		}
		auto score = stochastic.sScore(degree, edgesToCommunity[u], adjustedCommOut,
		                               adjustedExternalStubs);
		bool considerNode = isInCommunity[u] || score.first < scoreThreshold;
		if (considerNode)
			scores.emplace_back(score.first, score.second, u);
	}
	std::sort(scores.begin(), scores.end(), [](const ScoreStruct &a, const ScoreStruct &b) {
		return a.sScore < b.sScore;
	});
	return scores;
}

} /* namespace NetworKit */
