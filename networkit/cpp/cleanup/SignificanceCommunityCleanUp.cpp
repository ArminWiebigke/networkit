/*
 * SignificanceCommunityCleanUp.cpp
 *
 * Created: 2019-09-11
 * Author: Armin Wiebigke
 */

#include "SignificanceCommunityCleanUp.h"

namespace NetworKit {

using Community = SignificanceCommunityCleanUp::Community;

class SignificanceCommunityCleanUp::SingleCommunityCleanUp {
public:
	explicit SingleCommunityCleanUp(const Graph &graph, double scoreThreshold,
	                                double significanceThreshold);

	Community clean(const Community &inputCommunity);

	~SingleCommunityCleanUp() = default;

private:
	struct ScoreStruct {
		double sScore;
		double bootInterval;
		node candidate;

		ScoreStruct(double s, double b, node u) : sScore(s), bootInterval(b), candidate(u) {};
	};

	const Graph &graph;
	StochasticSignificance stochastic;
	// threshold to decide if a node is significant
	const double significanceThreshold;
	// threshold to discard candidates because they will most likely not be significant
	const double scoreThreshold;

	Community community;
	std::vector<count> kIn;         // kIn[u] == number of neighbors of node u that are in the community
	std::vector<int> isInCommunity; // isInCommunity[u] == is node u in the community
	std::vector<int> isCandidate;   // isCandidate[u] == is node u a candidate
	count outgoingCommunityStubs;   // number of edges that have one endpoint inside the community and one outside
	count totalCommunityStubs;      // the sum of the degrees of the nodes in the community
	count externalNodes;
	count externalStubs;

	std::vector<node> getCandidatesAndSetUpCalculation(bool includeNeighbors);

	void reset(const std::vector<node> &candidates);

	Community calculateSignificantNodes(const Community &inputCommunity, bool includeNeighbors);

	Community calculateSignificantCandidates(std::vector<ScoreStruct> scores) const;

	std::vector<SingleCommunityCleanUp::ScoreStruct>
	calculateCandidateScores(const std::vector<node> &candidates) const;

	void removeWorstNode(std::vector<ScoreStruct> scores);
};

SignificanceCommunityCleanUp::SignificanceCommunityCleanUp(const Graph &graph,
                                                           const Cover &cover,
                                                           double significanceThreshold,
                                                           double scoreThreshold,
                                                           double minOverlapRatio)
		: graph(graph),
		  cover(cover),
		  minOverlapRatio(minOverlapRatio) {
	singleCommunityCleanup.reset(
			new SingleCommunityCleanUp(graph, scoreThreshold, significanceThreshold));
}

SignificanceCommunityCleanUp::~SignificanceCommunityCleanUp() = default;

void SignificanceCommunityCleanUp::run() {
	hasRun = false;
	cleanedCommunities = Cover(graph.upperNodeIdBound());
	cleanAllCommunities();
	mergeDiscardedCommunities();
	hasRun = true;
}

Cover SignificanceCommunityCleanUp::getCover() {
	if (!hasFinished())
		throw std::runtime_error("Run the algorithm first!");
	return cleanedCommunities;
}

std::string SignificanceCommunityCleanUp::toString() const {
	return "SignificanceCommunityCleanUp";
}

bool SignificanceCommunityCleanUp::isParallel() const {
	return false;
}

void SignificanceCommunityCleanUp::cleanAllCommunities() {
	auto inputCommunities = cover.getSubsets();
	for (const Community &inputCommunity : inputCommunities) {
		auto cleanedCommunity = cleanCommunity(inputCommunity);
		if (cleanedCommunity.empty())
			discardedCommunities.push_back(inputCommunity);
		else
			cleanedCommunities.addSubset(cleanedCommunity);
	}
}

SignificanceCommunityCleanUp::Community
SignificanceCommunityCleanUp::cleanCommunity(const Community &inputCommunity) {
	Community cleanedCommunity = singleCommunityCleanup->clean(inputCommunity);
	bool changedDrastically = smallOverlap(inputCommunity, cleanedCommunity);
	if (changedDrastically)
		cleanedCommunity = {};
	return cleanedCommunity;
}

bool SignificanceCommunityCleanUp::smallOverlap(const Community &inputCommunity,
                                                const Community &cleanedCommunity) const {
	Community intersection;
	std::set_intersection(inputCommunity.begin(), inputCommunity.end(),
	                      cleanedCommunity.begin(), cleanedCommunity.end(),
	                      std::inserter(intersection, intersection.begin()));
	count largerSize = std::max(inputCommunity.size(), cleanedCommunity.size());
	double overlapRatio = double(intersection.size()) / largerSize;
	std::cout << inputCommunity.size() << " -> " << cleanedCommunity.size()
	          << (overlapRatio < minOverlapRatio ? " -> bad" : "") << std::endl;
	return overlapRatio < minOverlapRatio;
}

void SignificanceCommunityCleanUp::mergeDiscardedCommunities() {
	// TODO
}

Community
SignificanceCommunityCleanUp::SingleCommunityCleanUp::clean(const Community &inputCommunity) {
	Community firstPhaseResult = calculateSignificantNodes(inputCommunity, true);
	Community cleanedCommunity = calculateSignificantNodes(firstPhaseResult, false);
	return cleanedCommunity;
}

SignificanceCommunityCleanUp::SingleCommunityCleanUp::SingleCommunityCleanUp(
		const Graph &graph, double scoreThreshold, double significanceThreshold)
		: graph(graph),
		  scoreThreshold(scoreThreshold),
		  significanceThreshold(significanceThreshold),
		  kIn(graph.upperNodeIdBound()),
		  isInCommunity(graph.upperNodeIdBound()),
		  isCandidate(graph.upperNodeIdBound()),
		  stochastic(2 * graph.numberOfEdges() + graph.numberOfNodes()) {
}

Community
SignificanceCommunityCleanUp::SingleCommunityCleanUp::calculateSignificantNodes(
		const Community &inputCommunity, bool includeNeighbors) {
	community = inputCommunity;
	auto candidates = getCandidatesAndSetUpCalculation(includeNeighbors);
	auto testKIn = [&]() {
		for (node u : candidates) {
			assert(isInCommunity[u] == community.count(u));
		}
		std::map<node, count> kInTest;
		count commStubsTest = 0;
		count commOutTest = 0;
		for (node u : community) {
			graph.forNeighborsOf(u, [&](node v) {
				commStubsTest += 1;
				if (!isInCommunity[v])
					commOutTest += 1;
				if (isCandidate[v]) {
					kInTest[v] += 1;
				}
			});
		}
		for (node u : candidates) {
			assert(kIn[u] == kInTest[u]);
		}
		assert(outgoingCommunityStubs == commOutTest);
		assert(totalCommunityStubs == commStubsTest);
		return true;
	};
	assert(testKIn());
	Community cleanedCommunity;
	while (!community.empty()) {
		auto candidateScores = calculateCandidateScores(candidates);
		cleanedCommunity = calculateSignificantCandidates(candidateScores);
		bool communityIsSignificant = !cleanedCommunity.empty();
		if (communityIsSignificant)
			break;
		assert(testKIn());
		removeWorstNode(candidateScores);
		assert(testKIn());
	}
	reset(candidates);
	for (std::size_t i = 0; i < graph.upperNodeIdBound(); ++i) {
		assert(isInCommunity[i] == false);
		assert(isCandidate[i] == false);
		assert(kIn[i] == 0);

	}
	assert(community.empty());
	return cleanedCommunity;
}

std::vector<node>
SignificanceCommunityCleanUp::SingleCommunityCleanUp::getCandidatesAndSetUpCalculation(
		bool includeNeighbors) {
	std::vector<node> candidates;
	auto addCandidate = [&](node u) {
		candidates.push_back(u);
		isCandidate[u] = true;
	};
	for (node u : community) {
		isInCommunity[u] = true;
		addCandidate(u);
	}
	auto inOrigComm = isInCommunity;
	outgoingCommunityStubs = 0;
	totalCommunityStubs = 0;
	for (node u : community) {
		graph.forNeighborsOf(u, [&](node neighbor) {
			totalCommunityStubs += 1;
			if (!isInCommunity[neighbor])
				outgoingCommunityStubs += 1;
			if (isInCommunity[neighbor] || includeNeighbors) {
				if (!isCandidate[neighbor]) {
					assert(kIn[neighbor] == 0);
					addCandidate(neighbor);
				}
				kIn[neighbor] += 1;
			}
		});
	}
//	count externalNodes = graph.numberOfNodes() - community.size();
	externalNodes = graph.numberOfNodes();
	// TODO: Is the community part of external nodes? else we have more candidates than external nodes, which causes the order statistics to fail
	externalStubs = 2 * graph.numberOfEdges() - totalCommunityStubs;
	// TODO: external stubs with or without candidate degree
	return candidates;
}

void
SignificanceCommunityCleanUp::SingleCommunityCleanUp::reset(const std::vector<node> &candidates) {
	for (node u : candidates) {
		kIn[u] = 0;
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
void SignificanceCommunityCleanUp::SingleCommunityCleanUp::removeWorstNode(
		std::vector<ScoreStruct> scores) {
	auto worstNodeIt = std::find_if(scores.rbegin(), scores.rend(),
	                                [&](const ScoreStruct &s) { return isInCommunity[s.candidate]; });
	node worstNode = worstNodeIt->candidate;
	assert(isInCommunity[worstNode]);
	community.erase(worstNode);
	isInCommunity[worstNode] = false;
	// externalNodes += 1;
	// + kIn - kOut = + kIn - (k - kIn) = + 2 * kIn - k
	outgoingCommunityStubs += 2 * kIn[worstNode] - graph.degree(worstNode);
	totalCommunityStubs -= graph.degree(worstNode);
	graph.forNeighborsOf(worstNode, [&](node u) {
		if (isCandidate[u]) {
			assert(kIn[u] > 0);
			kIn[u] -= 1;
		}
	});
}

Community
SignificanceCommunityCleanUp::SingleCommunityCleanUp::calculateSignificantCandidates(
		std::vector<ScoreStruct> scores) const {
	int position = 1;
	int significantNodesCount = 0;
	for (auto scoreStruct : scores) {
		double score = scoreStruct.sScore;
		// significance is the probability that Omega_{position} <= score == Phi(significance)
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

std::vector<SignificanceCommunityCleanUp::SingleCommunityCleanUp::ScoreStruct>
SignificanceCommunityCleanUp::SingleCommunityCleanUp::calculateCandidateScores(
		const std::vector<node> &candidates) const {
	std::vector<ScoreStruct> scores;
	for (node u : candidates) {
		// TODO: How to handle external stubs if candidate is/is not currently in community
		count adjustedExternalStubs = externalStubs;
		// Calculate scores as if the node is not part of the community
		count adjustedCommOut = outgoingCommunityStubs;
		if (isInCommunity[u])
			adjustedCommOut += 2 * kIn[u] - graph.degree(u);
		auto score = stochastic.sScore(graph.degree(u), kIn[u], adjustedCommOut,
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
