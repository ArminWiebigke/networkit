/*
 * SignificanceCommunityCleanUp.cpp
 *
 * Created: 2019-09-11
 * Author: Armin Wiebigke
 */

#include "SignificanceCommunityCleanUp.h"

namespace NetworKit {

SignificanceCommunityCleanUp::SignificanceCommunityCleanUp(const Graph &graph, const Cover &cover)
		: graph(graph),
		  cover(cover),
		  gamma(0.5),
		  singleCommunityCleanup(graph) {

}

void SignificanceCommunityCleanUp::run() {
	init();
	checkComms();
	mergeDiscarded();
	hasRun = true;
}

void SignificanceCommunityCleanUp::init() {
	hasRun = false;
	resultCover = Cover(graph.upperNodeIdBound());
}

Cover SignificanceCommunityCleanUp::getCover() {
	if (!hasFinished())
		throw std::runtime_error("Run the algorithm first!");
	return resultCover;
}

std::string SignificanceCommunityCleanUp::toString() const {
	return "SignificanceCommunityCleanUp";
}

bool SignificanceCommunityCleanUp::isParallel() const {
	return false;
}

void SignificanceCommunityCleanUp::checkComms() {
	auto communities = cover.getSubsets();
	for (const auto &community : communities) {
		auto cleanedComm = cleanCommunity(community);
		if (!cleanedComm.empty())
			resultCover.addSubset(cleanedComm);
		else
			discardedComms.push_back(community);
	}
}

SignificanceCommunityCleanUp::Community
SignificanceCommunityCleanUp::cleanCommunity(const Community &community) {
	Community cleanedComm = singleCommunityCleanup.clean(community);

	// Make sure the community did not change drastically
	Community intersect;
	std::set_intersection(community.begin(), community.end(),
	                      cleanedComm.begin(), cleanedComm.end(),
	                      std::inserter(intersect, intersect.begin()));
	double overlapRatio = double(intersect.size()) / std::max(community.size(), cleanedComm.size());
	std::cout << community.size() << " -> " << cleanedComm.size()
	          << (overlapRatio < gamma ? " -> bad" : "") << std::endl;
	if (overlapRatio < gamma)
		cleanedComm = {};

	return cleanedComm;
}

void SignificanceCommunityCleanUp::mergeDiscarded() {
	// TODO
}

SignificanceCommunityCleanUp::SingleCommunityCleanup::SingleCommunityCleanup(const Graph &graph)
		: graph(graph),
		  kIn(graph.upperNodeIdBound()),
		  inComm(graph.upperNodeIdBound()),
		  stochastic(2 * graph.numberOfEdges() + graph.numberOfNodes()) {

}

SignificanceCommunityCleanUp::Community SignificanceCommunityCleanUp::SingleCommunityCleanup::clean(
		const Community &community) {
	Community cleanedComm = community;
	cleanedComm = checkCommSignificance(community, true);
	cleanedComm = checkCommSignificance(cleanedComm, false);
	return cleanedComm;
}

SignificanceCommunityCleanUp::Community
SignificanceCommunityCleanUp::SingleCommunityCleanup::checkCommSignificance(
		Community community, bool checkNeighbors) {
	Community curCommunity;
	Community originalComm = community;
	Community cleanedComm = {};
	std::vector<node> candidates;
	for (node u : community) {
		candidates.push_back(u);
		inComm[u] = true;
	}
	auto inOrigComm = inComm;
	count commOut = 0; // stubs outgoing from the community
	count commStubs = 0; // the sum of the degrees of the nodes in the community
	for (node u : community) {
		graph.forNeighborsOf(u, [&](node v) {
			commStubs += 1;
			if (!inComm[v])
				commOut += 1;
			if (inComm[v] || checkNeighbors) {
				if (kIn[v] == 0 && !inComm[v])
					candidates.push_back(v);
				kIn[v] += 1;
			}
		});
	}
//	count externalNodes = graph.numberOfNodes() - community.size();
	count externalNodes = graph.numberOfNodes();
	// TODO: Is the community part of external nodes? else we have more candidates than external nodes, which causes the order statistics to fail
	count externalStubs = 2 * graph.numberOfEdges() - commStubs;
	// TODO: external stubs with or without candidate degree

	auto testKIn = [&]() {
		for (node u : candidates) {
			assert(inComm[u] == community.count(u));
		}
		std::map<node, count> kInTest;
		count commStubsTest = 0;
		count commOutTest = 0;
		for (node u : community) {
			graph.forNeighborsOf(u, [&](node v) {
				commStubsTest += 1;
				if (!inComm[v])
					commOutTest += 1;
				if (inOrigComm[v] || checkNeighbors) {
					kInTest[v] += 1;
				}
			});
		}
		for (node u : candidates) {
			assert(kIn[u] == kInTest[u]);
		}
		assert(commOut == commOutTest);
		assert(commStubs == commStubsTest);
		return true;
	};
	assert(testKIn());

	while (!community.empty()) {
		// Calculate s-Scores
		std::vector<ScoreStruct> scores;
		for (node u : candidates) {
			// TODO: How to handle external stubs if candidate is/is not currently in community
			count tmpExternalStubs = externalStubs;
			count tmpCommOut = inComm[u] ? commOut + 2 * kIn[u] - graph.degree(u) : commOut;
			auto score = stochastic.sScore(graph.degree(u), kIn[u], tmpCommOut, tmpExternalStubs);
			if (score.first < scoreThreshold ||
			    inComm[u]) // only consider neighbors with a good score
				scores.emplace_back(score.first, score.second, u);
		}
		std::sort(scores.begin(), scores.end(), [](ScoreStruct a, ScoreStruct b) {
			return a.sScore < b.sScore;
		});

		// Find significant candidates
		cleanedComm = significantCandidates(scores, externalNodes);

		// If the group is not empty, the module is considered significant
		if (!cleanedComm.empty()) {
			break;
		}

		assert(testKIn());

		// Remove worst node from community
		auto it = std::find_if(scores.rbegin(), scores.rend(),
		                       [&](const ScoreStruct &s) { return inComm[s.u]; });
		node worstNode = it->u;
		assert(inComm[worstNode]);
		community.erase(worstNode);
		inComm[worstNode] = false;
//		externalNodes += 1;
		// + kIn - kOut = + kIn - (k - kIn) = + 2 * kIn - k
		commOut += 2 * kIn[worstNode] - graph.degree(worstNode);
		commStubs -= graph.degree(worstNode);
		graph.forNeighborsOf(worstNode, [&](node v) {
			if (inOrigComm[v] || checkNeighbors) {
				assert(kIn[v] > 0);
				kIn[v] -= 1;
			}
		});

		assert(testKIn());
	}

	// Reset data structures
	for (node u : candidates)
		kIn[u] = 0;
	for (node u : community)
		inComm[u] = false;
	return cleanedComm;
}

SignificanceCommunityCleanUp::Community
SignificanceCommunityCleanUp::SingleCommunityCleanup::significantCandidates(
		const std::vector<ScoreStruct> &scores, count externalNodes) {
	int pos = 1; // position for the order statistic
	int sigCnt = 0; // sigCnt tells how many nodes should be included into the cluster

	for (auto scoreStruct : scores) {
		double score = scoreStruct.sScore;
		/*
		 significance is the probability that Omega_{pos} <= score == Phi(significance)
		 If this is low, the score of the node is better than expected in the null model,
		 so we can assume that the node is part of the community
		*/
		double significance = stochastic.orderStatistic(score, externalNodes, pos);
		if (significance < sigThreshold) {
			sigCnt = pos;
		} else {
			if (sigCnt != 0) {
				// Stop after finding significant nodes and then not significant ones
				break;
			}
		}
		++pos;
	}

	// Insert nodes in the cleaned group
	Community significantNodes;
	for (count i = 0; i < sigCnt; ++i) {
		significantNodes.insert(scores[i].u);
	}
	return significantNodes;
}

} /* namespace NetworKit */
