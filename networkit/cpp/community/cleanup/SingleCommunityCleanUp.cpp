/*
 * SingleCommunityCleanup.cpp
 *
 * Created: 2019-09-26
 * Author: Armin Wiebigke
 */
#include "SingleCommunityCleanUp.h"

namespace NetworKit {

using Community = SingleCommunityCleanUp::Community;

SingleCommunityCleanUp::SingleCommunityCleanUp(const Graph &graph,
                                               const StochasticDistribution &stochasticDistribution,
                                               double scoreThreshold,
                                               double significanceThreshold, double minOverlapRatio)
		: graph(graph),
		  significanceThreshold(significanceThreshold),
		  scoreThreshold(scoreThreshold),
		  minOverlapRatio(minOverlapRatio),
		  edgesToCommunity(graph.upperNodeIdBound()),
		  isInCommunity(graph.upperNodeIdBound()),
		  isInOriginalCommunity(graph.upperNodeIdBound()),
		  isCandidate(graph.upperNodeIdBound()),
		  significanceCalculator({stochasticDistribution}) {
}

Community
SingleCommunityCleanUp::clean(const Community &inputCommunity) {
	Community firstPhaseResult = firstPhase(inputCommunity);
	Community cleanedCommunity = secondPhase(firstPhaseResult);
	bool changedDrastically = smallOverlap(inputCommunity, cleanedCommunity);
	if (changedDrastically)
		cleanedCommunity = {};
	return cleanedCommunity;
}

Community
SingleCommunityCleanUp::secondPhase(const Community &firstPhaseResult) {
	return calculateSignificantNodes(firstPhaseResult, true);
}

Community
SingleCommunityCleanUp::firstPhase(const Community &inputCommunity) {
	return calculateSignificantNodes(inputCommunity, false);
}

/**
 * Returns a community of significant nodes
 * @param onlyUseOriginalCommunity If false, consider neighbors of the input community and the nodes
 * of the input community. If true, consider only the nodes of the input community.
 */
Community
SingleCommunityCleanUp::calculateSignificantNodes(
		const Community &inputCommunity, bool onlyUseOriginalCommunity) {
	community = inputCommunity;
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
	for (node u : community) {
		isInCommunity.insert(u, true);
		isInOriginalCommunity.insert(u, true);
	}
	auto tryAddCandidate = [&](node u) {
		if (!isInCommunity[u] && !isCandidate[u]) {
			isCandidate.insert(u, true);
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
				if (edgesToCommunity[neighbor] == 0)
					edgesToCommunity.insert(neighbor, 0);
				edgesToCommunity[neighbor] += 1;
			}
		});
	}

	externalNodes = graph.numberOfNodes() - community.size();
	externalStubs = 2 * graph.numberOfEdges() - totalCommunityStubs;
}

// remove the node with the worst (= highest) score from the community
void SingleCommunityCleanUp::removeWorstNode(
		std::vector<ScoreStruct> &internalScores) {
	assert(community.size() > 0);
	assert(internalScores.size() > 0);

	auto removeNode = [&](node nodeToRemove) {
		assert(isInCommunity[nodeToRemove]);
		community.erase(nodeToRemove);
		isInCommunity[nodeToRemove] = false;
		candidates.push_back(nodeToRemove);
		isCandidate.insert(nodeToRemove, true);
		externalNodes += 1;
		assert(externalNodes == graph.numberOfNodes() - community.size());
		count degree = graph.degree(nodeToRemove);
		outgoingCommunityStubs += 2 * edgesToCommunity[nodeToRemove] - degree;
		totalCommunityStubs -= degree;
		externalStubs += degree;
		graph.forNeighborsOf(nodeToRemove, [&](node u) {
			if (isCandidate[u] || isInCommunity[u]) {
				assert(edgesToCommunity[u] > 0);
				edgesToCommunity[u] -= 1;
			}
		});
	};

	if (internalScores.size() < 20) {
		// TODO: Calculate the score of the nodes inside the community separately
		auto worstNodeIt = std::max_element(internalScores.begin(), internalScores.end(),
						    [](const ScoreStruct &a, const ScoreStruct &b) {
							    return a.rScore < b.rScore;
						    });
		node worstNode = worstNodeIt->nodeId;
		removeNode(worstNode);
	} else {
		count nodesToRemove = internalScores.size() / 10;
		std::partial_sort(internalScores.begin(), internalScores.begin() + nodesToRemove, internalScores.end(),
						    [](const ScoreStruct &a, const ScoreStruct &b) {
							    return a.rScore > b.rScore;
						    });
		for (count i = 0; i < nodesToRemove; ++i) {
			removeNode(internalScores[i].nodeId);
		}
	}
}

std::vector<SingleCommunityCleanUp::ScoreStruct>
SingleCommunityCleanUp::calculateInternalScores() {
	std::vector<ScoreStruct> internalScores;
	internalScores.reserve(community.size());
	for (node u : community) {
		count degree = graph.degree(u);
		// Calculate s-Score as if the node was not part of the community
		count adjustedOutgoingStubs = outgoingCommunityStubs + 2 * edgesToCommunity[u] - degree;
		count adjustedExternalStubs = externalStubs + degree;
		double score = significanceCalculator.rScore(
				degree, edgesToCommunity[u], adjustedOutgoingStubs, adjustedExternalStubs);
		internalScores.emplace_back(score, u);
	}
	return internalScores;
}


std::vector<SingleCommunityCleanUp::ScoreStruct>
SingleCommunityCleanUp::calculateCandidateScores() {
	std::vector<ScoreStruct> candidateScores;
	for (node u : candidates) {
		assert(!isInCommunity[u]);
		double score = significanceCalculator.rScore(graph.degree(u), edgesToCommunity[u],
		                                 outgoingCommunityStubs, externalStubs);
		if (score < scoreThreshold)
			candidateScores.emplace_back(score, u);
	}
	std::sort(candidateScores.begin(), candidateScores.end(),
	          [](const ScoreStruct &a, const ScoreStruct &b) {
		          return a.rScore < b.rScore;
	          });
	return candidateScores;
}

namespace {
// TODO: What is the reason behind this?
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
}

std::vector<node>
SingleCommunityCleanUp::findSignificantCandidates(const std::vector<ScoreStruct> &scores) {
	int position = 1;
	int significantNodesCount = 0;
	double threshold = significanceThreshold / fitted_exponent(externalNodes);
	for (auto scoreStruct : scores) {
		double score = scoreStruct.rScore;
		// significance is the probability Omega_{position}(score, externalNodes)
		double significance = significanceCalculator.orderStatistic(score, externalNodes, position);
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
	significantNodes.reserve(significantNodesCount);
	for (index i = 0; i < significantNodesCount; ++i) {
		significantNodes.push_back(scores[i].nodeId);
	}
	return significantNodes;
}

void
SingleCommunityCleanUp::reset() {
	isCandidate.reset();
	edgesToCommunity.reset();
	isInCommunity.reset();
	isInOriginalCommunity.reset();
#ifndef NDEBUG
	for (index i = 0; i < graph.upperNodeIdBound(); ++i) {
		assert(edgesToCommunity[i] == 0);
		assert(isCandidate[i] == 0);
		assert(isInCommunity[i] == 0);
		assert(isInOriginalCommunity[i] == 0);
		assert(edgesToCommunity[i] == 0);
	}
#endif
	community.clear();
	candidates.clear();
}

bool SingleCommunityCleanUp::smallOverlap(const Community &inputCommunity,
                                          const Community &cleanedCommunity) const {
	std::vector<node> overlap;
	std::set_intersection(inputCommunity.begin(), inputCommunity.end(),
	                      cleanedCommunity.begin(), cleanedCommunity.end(),
	                      std::back_inserter(overlap));
	count largerSize = std::max(inputCommunity.size(), cleanedCommunity.size());
	double overlapRatio = (double) overlap.size() / largerSize;
//	std::cout << inputCommunity.size() << " -> " << cleanedCommunity.size()
//	          << (overlapRatio < minOverlapRatio ? " -> bad" : "") << std::endl;
	return overlapRatio < minOverlapRatio;
}

} /* namespace NetworKit */
