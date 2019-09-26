/*
 * MergeCommunities.h
 *
 * Created: 2019-09-26
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_MERGECOMMUNITIES_H
#define NETWORKIT_MERGECOMMUNITIES_H

#include "../graph/Graph.h"
#include "SingleCommunityCleanUp.h"
#include "../structures/Partition.h"

namespace NetworKit {

class MergeCommunities {
public:
	using Community = std::set<node>;

	MergeCommunities(const Graph &graph, std::set<Community> discardedCommunities,
	                 SingleCommunityCleanUp &singleCommunityCleanUp)
			: graph(graph),
			  discardedCommunities(std::move(discardedCommunities)),
			  stochastic(graph.numberOfNodes() + 2 * graph.numberOfEdges()),
			  singleCommunityCleanUp(singleCommunityCleanUp) {
	}

	void run();

	std::vector<Community> & getCleanedCommunities();

private:
	const Graph &graph;
	std::vector<Community> cleanedCommunities;
	std::set<Community> discardedCommunities;
	Graph discardedCommunitiesGraph;
	Partition mergedCommunities;
	std::vector<std::set<node>> coarseToFineMapping;
	StochasticSignificance stochastic;
	std::vector<count> cOut;
	std::vector<count> cTotal;
	count totalStubs;
	SingleCommunityCleanUp &singleCommunityCleanUp;

	void createDiscardedCommunitiesGraph();

	void localOptimization();

	void checkMergedCommunities();

	bool tryLocalMove(node u);
};


} /* namespace NetworKit */


#endif //NETWORKIT_MERGECOMMUNITIES_H
