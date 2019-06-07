/*
 * SignificanceCleanup.cpp
 *
 * Created: 2019-05-02
 * Author: Armin Wiebigke
 */

#include "SignificanceCleanup.h"

namespace NetworKit {

SignificanceCleanup::SignificanceCleanup(const Graph &G, const Cover &input)
		: G(G), cover(input), Algorithm() {
}

// TODO: Use this?
void SignificanceCleanup::run() {
	std::vector<std::set<node>> communities(cover.upperBound());
	G.forNodes([&](node u){
		for (index c : cover.subsetsOf(u)) {
			communities[c].insert(u);
		}
	});

	for (auto &c : communities) {
		cleanupCommunity(c);
	}

	cleanedCover = Cover(G.upperNodeIdBound());
	cleanedCover.setUpperBound(communities.size());
	for (size_t c_id = 0; c_id < communities.size(); ++c_id) {
		auto c = communities[c_id];
		for (node u : c) {
			cleanedCover.addToSubset(c_id, u);
		}
	}

	hasRun = true;
}

Cover SignificanceCleanup::getCover() {
	return cleanedCover;
}

void SignificanceCleanup::cleanupCommunity(std::set<node> &community) {
	while (!community.empty()) {
		auto rScores = calculateRScores(community);

		// Get worst node in the community
		auto maxScore = getMaxRScore(rScores);
		node worstNode = maxScore.first;

		// Abort the cleanup for this community if the worst node is significant,
		// else remove the node and continue the cleanup
		if (isSignificant(worstNode))
			break;
		else
			community.erase(worstNode);
	}
}

std::map<node, double> SignificanceCleanup::calculateRScores(
		std::set<node> &community) {
	std::map<node, double> rScores;

	for(node u : community) {
		rScores[u] = calculateRScore(u);
	}

	return rScores;
}

std::pair<node, double> SignificanceCleanup::getMaxRScore(
		std::map<node, double> &significances) {
	node maxNode = 0;
	double maxScore = -1.0;
	for (auto &pair : significances) {
		double score = pair.second;
		if (score > maxScore) {
			maxScore = score;
			maxNode = pair.first;
		}
	}
	return {maxNode, maxScore};
}

double SignificanceCleanup::calculateRScore(node u) {
	return 0;
}

bool SignificanceCleanup::isSignificant(node u) {
	return false;
}

} /* namespace NetworKit */