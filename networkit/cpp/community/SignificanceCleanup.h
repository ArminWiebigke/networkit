/*
 * SignificanceCleanup.h
 *
 * Created: 2019-05-02
 * Author: Armin Wiebigke
 */

#ifndef SIGNIFICANCECLEANUP_H
#define SIGNIFICANCECLEANUP_H

#include "../base/Algorithm.h"
#include "../graph/Graph.h"
#include "../structures/Cover.h"

namespace NetworKit {

class SignificanceCleanup : public Algorithm {
public:
	SignificanceCleanup(const Graph &G, const Cover &input);

	void run() override;

	Cover getCover();

private:
	const Graph &G;
	Cover cover, cleanedCover;

	void cleanupCommunity(std::set<node> &community);

	std::map<node, double> calculateRScores(std::set<node> &community);

	std::pair<node, double> getMaxRScore(std::map<node, double> &significances);

	double calculateRScore(node u);

	bool isSignificant(node u);

};

} /* namespace NetworKit */

#endif //SIGNIFICANCECLEANUP_H
