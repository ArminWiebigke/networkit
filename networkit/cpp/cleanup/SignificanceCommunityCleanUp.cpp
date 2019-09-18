/*
 * SignificanceCommunityCleanUp.cpp
 *
 * Created: 2019-09-11
 * Author: Armin Wiebigke
 */

#include "SignificanceCommunityCleanUp.h"

namespace NetworKit {

SignificanceCommunityCleanUp::SignificanceCommunityCleanUp(const Graph &graph, const Cover &cover)
		: graph(graph), cover(cover), resultCover(graph.upperNodeIdBound()), gamma(0.5),
		  kIn(graph.upperNodeIdBound(), 0), inComm(graph.upperNodeIdBound()),
		  stochastic(2 * graph.numberOfEdges() + graph.numberOfNodes()) {

}

void SignificanceCommunityCleanUp::run() {
	checkComms();
	mergeDiscarded();

	hasRun = true;
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

Cover SignificanceCommunityCleanUp::getCover() {
	return resultCover;
}

std::string SignificanceCommunityCleanUp::toString() const {
	return "SignificanceCommunityCleanUp";
}

bool SignificanceCommunityCleanUp::isParallel() const {
	return false;
}

SignificanceCommunityCleanUp::Community
SignificanceCommunityCleanUp::cleanCommunity(const Community &community) {
	if (community.size() < 2)
		return {};

	Community cleanedComm;
	cleanedComm = checkCommSignificance(community, true);
	if (cleanedComm.empty())
		return {};
	cleanedComm = checkCommSignificance(cleanedComm, false);
	if (cleanedComm.empty())
		return {};

	// Make sure the community did not change drastically
	std::set<int> intersect;
	std::set_intersection(community.begin(), community.end(),
	                      cleanedComm.begin(), cleanedComm.end(),
	                      std::inserter(intersect, intersect.begin()));
	double overlapRatio = double(intersect.size()) / std::max(community.size(), cleanedComm.size());
	if (overlapRatio < gamma)
		return {};

	return cleanedComm;
}

void SignificanceCommunityCleanUp::mergeDiscarded() {
	// TODO
}

SignificanceCommunityCleanUp::Community
SignificanceCommunityCleanUp::checkCommSignificance(Community community, bool checkNeighbors) {
	Community originalComm = community;
	Community cleanedComm = {};
	std::vector<node> candidates;
	for (std::size_t i = 0; i < inComm.size(); ++i) {
	    assert(inComm[i] == false);
	}
	for (node u : community)
		inComm[u] = true;
	count commOut = 0;
	count commStubs = 0;
	for (std::size_t i = 0; i < kIn.size(); ++i) {
	    assert(kIn[i] == 0);
	}
	for (node u : community) {
		graph.forNeighborsOf(u, [&](node v) {
			commStubs += 1;
			if (!inComm[v])
				commOut += 1;
			if (inComm[v] || checkNeighbors) {
				if (kIn[v] == 0)
					candidates.push_back(v);
				kIn[v] += 1;
			}
		});
	}
	auto ogKIn = kIn;
//	count externalNodes = graph.numberOfNodes() - community.size();
	count externalNodes = graph.numberOfNodes(); // TODO: community part of external nodes? else
	//   we have more candidates than external nodes, which causes the order statistics to fail
	count externalStubs = 2 * graph.numberOfEdges() - commStubs;
	// TODO: external stubs with or without candidate degree

	while (true) {
		// Calculate s-Scores
		std::vector<ScoreStruct> scores;
		for (node u : candidates) {
			// TODO: How to handle external stubs if candidate is/is not currently in community
			count tmpExternalStubs = externalStubs;
			count tmpCommOut = inComm[u] ? commOut + 2 * kIn[u] - graph.degree(u) : commOut;
			auto score = stochastic.sScore(kIn[u], tmpCommOut, tmpExternalStubs, graph.degree(u));
			if (score.first < 0.99) // only consider candidates with a good score
				scores.emplace_back(score.first, score.second, u);
		}
		std::sort(scores.begin(), scores.end(), [](ScoreStruct a, ScoreStruct b) {
			return a.sScore < b.sScore;
		});
		if (scores.empty()) {
			cleanedComm = {};
			break;
		}

		// TODO: Find significant candidates
		// cleanedComm =
		cleanedComm = significantCandidates(scores, externalNodes);

		// If the group is not empty, the module is considered significant
		if (!cleanedComm.empty()) {
			break;
		}

		// Remove worst node from community
		node worstNode = scores.back().u;
		community.erase(worstNode);
		inComm[worstNode] = false;
		externalNodes += 1;
		// + kIn - kOut = + kIn - (k - kIn) = + 2 * kIn - k
		commOut +=  2 * kIn[worstNode] - graph.degree(worstNode);
		graph.forNeighborsOf(worstNode, [&](node v) {
			if (inComm[v] || checkNeighbors) {
				assert(kIn[v] > 0);
				kIn[v] -= 1;
			}
		});

		for (auto b : kIn)
			assert(b < 1e6);

		// If all nodes have been removed, the module is insignificant
		if (community.empty()) {
			cleanedComm = {};
			break;
		}
	}

	// Reset data structures
	for (node u : candidates)
		kIn[u] = 0;
	for (node u : candidates)
		assert(kIn[u] == 0);
	for (auto b : kIn)
		assert(b == 0);
	for (node u : community)
		inComm[u] = false;

	return cleanedComm;
}

SignificanceCommunityCleanUp::Community
SignificanceCommunityCleanUp::significantCandidates(const std::vector<ScoreStruct> &scores,
                                                    count externalNodes) {
	int pos = 1;
	int sigCnt = 0; // sigCnt tells how many nodes should be included into the cluster

	for (auto scoreStruct : scores) {
		double score = scoreStruct.sScore;
		/*
		 c_pos is the probability that Omega_{pos} <= score == Phi(c_pos)
		 If this is low, the score of the node is better than expected in the null model,
		 so we can assume that the node is part of the community
		*/
		//TODO: Is this really the correct value
		double c_pos = stochastic.orderStatistic(score, externalNodes, pos);

		if (c_pos < sigThreshold) {
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
//	Community significantNodes(scores.begin(), scores.begin() + sigCnt);
	return significantNodes;
}

} /* namespace NetworKit */
