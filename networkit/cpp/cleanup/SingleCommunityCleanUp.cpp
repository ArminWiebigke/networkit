/*
 * SingleCommunityCleanup.cpp
 *
 * Created: 2019-09-26
 * Author: Armin Wiebigke
 */
#include "SingleCommunityCleanUp.h"

namespace NetworKit {

using Community = SingleCommunityCleanUp::Community;

SingleCommunityCleanUp::SingleCommunityCleanUp(const Graph &graph, double scoreThreshold,
                                               double significanceThreshold, double minOverlapRatio)
		: graph(graph),
		  scoreThreshold(scoreThreshold),
		  significanceThreshold(significanceThreshold),
		  minOverlapRatio(minOverlapRatio),
		  edgesToCommunity(graph.upperNodeIdBound()),
		  isInCommunity(graph.upperNodeIdBound()),
		  isInOriginalCommunity(graph.upperNodeIdBound()),
		  isCandidate(graph.upperNodeIdBound()),
		  stochastic(2 * graph.numberOfEdges() + graph.numberOfNodes()) {
}

Community
SingleCommunityCleanUp::clean(const Community &inputCommunity) {
	Community firstPhaseResult = calculateSignificantNodes(inputCommunity, false);
	Community cleanedCommunity = calculateSignificantNodes(firstPhaseResult, true);
	bool changedDrastically = smallOverlap(inputCommunity, cleanedCommunity);
	if (changedDrastically)
		cleanedCommunity = {};
	return cleanedCommunity;
}

Community
SingleCommunityCleanUp::calculateSignificantNodes(
		const Community &inputCommunity, bool onlyUseOriginalCommunity) {
	originalCommunity = inputCommunity;
	community = originalCommunity;
	getCandidatesAndSetUpCalculation(onlyUseOriginalCommunity);
	std::vector<node> significantNeighbors;
	while (!community.empty()) {
		auto internalScores = calculateInternalScores();
		auto candidateScores = calculateCandidateScores();
		significantNeighbors = findSignificantCandidates(candidateScores);
		bool communityIsSignificant = !significantNeighbors.empty();
		if (communityIsSignificant)
			break;
		removeWorstNode(internalScores);
	}
	Community significantNodes = community;
	for (node u : significantNeighbors) {
		significantNodes.insert(u);
	}
	reset();
	return significantNodes;
}

void
SingleCommunityCleanUp::getCandidatesAndSetUpCalculation(bool onlyUseOriginalCommunity) {
	candidates.clear();
	for (node u : community) {
		isInCommunity[u] = true;
		isInOriginalCommunity[u] = true;
		assert(originalCommunity.count(u) == 1);
	}
	auto tryAddCandidate = [&](node u) {
		if (!isInCommunity[u] && !isCandidate[u]) {
			isCandidate[u] = true;
			candidates.push_back(u);
			assert(edgesToCommunity[u] == 0);
		}
	};
	outgoingCommunityStubs = 0;
	totalCommunityStubs = 0;
	for (node u : community) {
		graph.forNeighborsOf(u, [&](node neighbor) {
			totalCommunityStubs += 1;
			if (!isInCommunity[neighbor])
				outgoingCommunityStubs += 1;
			bool validCandidate = isInOriginalCommunity[neighbor] || !onlyUseOriginalCommunity;
			if (validCandidate) {
				tryAddCandidate(neighbor);
				edgesToCommunity[neighbor] += 1;
			}
		});
	}

	externalNodes = graph.numberOfNodes() - community.size();
	externalStubs = 2 * graph.numberOfEdges() - totalCommunityStubs;
}

void
SingleCommunityCleanUp::reset() {
	for (node u : candidates) {
		edgesToCommunity[u] = 0;
		isCandidate[u] = 0;
	}
	for (node u : community) {
		edgesToCommunity[u] = 0;
		isInCommunity[u] = false;
	}
	for (node u : originalCommunity) {
		isInOriginalCommunity[u] = false;
		assert(edgesToCommunity[u] == 0);
	}
	for (index i = 0; i < graph.upperNodeIdBound(); ++i) {
		assert(edgesToCommunity[i] == 0);
		assert(isCandidate[i] == 0);
		assert(isInCommunity[i] == 0);
		assert(isInOriginalCommunity[i] == 0);
		assert(edgesToCommunity[i] == 0);
	}
}

// remove the node with the worst (= highest) score from the community
void SingleCommunityCleanUp::removeWorstNode(
		std::vector<ScoreStruct> internalScores) {
	assert(community.size() > 0);
	assert(internalScores.size() > 0);
	// TODO: Calculate the score of the nodes inside the community separately
	auto worstNodeIt = std::max_element(internalScores.begin(), internalScores.end(),
	                                    [](const ScoreStruct &a, const ScoreStruct &b) {
		                                    return a.sScore < b.sScore;
	                                    });
	node worstNode = worstNodeIt->nodeId;
	assert(isInCommunity[worstNode]);
	community.erase(worstNode);
	isInCommunity[worstNode] = false;
	candidates.push_back(worstNode);
	isCandidate[worstNode] = true;
	externalNodes += 1;
	assert(externalNodes == graph.numberOfNodes() - community.size());
	count degree = graph.degree(worstNode);
	outgoingCommunityStubs += 2 * edgesToCommunity[worstNode] - degree;
	totalCommunityStubs -= degree;
	externalStubs += degree;
	graph.forNeighborsOf(worstNode, [&](node u) {
		if (isCandidate[u] || isInCommunity[u]) {
			assert(edgesToCommunity[u] > 0);
			edgesToCommunity[u] -= 1;
		}
	});
}

double fitted_exponent(int N) {
	double l = log(double(N));

	if (N > 100)
		return 4.2 * l - 8.5;

	if (N > 30)
		return 3.5 * l - 5.5;

	if (N > 7)
		return 2.5 * l - 2;

	if (N > 1)
		return 1.3 * l + 0.1;

	return 1;
}

std::vector<node>
SingleCommunityCleanUp::findSignificantCandidates(std::vector<ScoreStruct> scores) const {
	int position = 1;
	int significantNodesCount = 0;
	for (auto scoreStruct : scores) {
		double score = scoreStruct.sScore;
		// significance is the probability Omega_{position}(score, externalNodes)
		double significance = stochastic.orderStatistic(score, externalNodes, position);
//		double threshold = significanceThreshold;
		double threshold = significanceThreshold / fitted_exponent(externalNodes);
		if (significance < threshold) {
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

	std::vector<node> significantNodes;
	for (index i = 0; i < significantNodesCount; ++i) {
		significantNodes.push_back(scores[i].nodeId);
	}
	return significantNodes;
}


std::vector<SingleCommunityCleanUp::ScoreStruct>
SingleCommunityCleanUp::calculateCandidateScores() const {
	std::vector<ScoreStruct> scores;
	for (node u : candidates) {
		assert(!isInCommunity[u]);
		auto score = stochastic.sScore(graph.degree(u), edgesToCommunity[u], outgoingCommunityStubs,
		                               externalStubs);
		if (score.first < scoreThreshold)
			scores.emplace_back(score.first, score.second, u);
	}
	std::sort(scores.begin(), scores.end(), [](const ScoreStruct &a, const ScoreStruct &b) {
		return a.sScore < b.sScore;
	});
	return scores;
}

std::vector<SingleCommunityCleanUp::ScoreStruct>
SingleCommunityCleanUp::calculateInternalScores() {
	std::vector<ScoreStruct> internalScores;
	for (node u : community) {
		double score, bootInterval;
		count degree = graph.degree(u);
		// Calculate s-Score as if the node was not part of the community
		count adjustedOutgoingStubs = outgoingCommunityStubs + 2 * edgesToCommunity[u] - degree;
		count adjustedExternalStubs = externalStubs + degree;
		std::tie(score, bootInterval) = stochastic.sScore(
				degree, edgesToCommunity[u], adjustedOutgoingStubs, adjustedExternalStubs);
		internalScores.emplace_back(score, bootInterval, u);
	}
	return internalScores;
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

} /* namespace NetworKit */
