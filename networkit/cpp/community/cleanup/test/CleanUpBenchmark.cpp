/*
 * CleanUpBenchmark.cpp
 *
 * Created: 2019-09-29
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>

#include "../../../auxiliary/Timer.h"
#include "../../../graph/Graph.h"
#include "../../../io/EdgeListReader.h"
#include "../../egosplitting/EgoSplitting.h"
#include "../SignificanceCommunityCleanUp.h"

namespace NetworKit {

class CleanUpBenchmark : public testing::Test {
};

TEST_F(CleanUpBenchmark, benchCommunityCleanup) {
//	METISGraphReader graphReader;
//	Graph G = graphReader.read("../input/10_clusters.graph");
//                             G);
	EdgeListReader reader(' ', 0);
	Graph G = reader.read("../input/lfr_om3.graph");
//		Graph G = reader.read("/home/armin/graphs/email-Eu-core.txt");
//		G.removeSelfLoops();

	Aux::Timer timer;
	timer.start();
	EgoSplitting algo(G);
	algo.run();
	Cover cover = algo.getCover();
	timer.stop();
	std::cout << "egosplitting took " << timer.elapsedMilliseconds() << "ms" << std::endl;

	timer.start();
	SignificanceCommunityCleanUp cleanUp(G, cover, 0.1, 0.1, 0.5);
	cleanUp.run();
	Cover cleanedCover = cleanUp.getCover();
	timer.stop();
	std::cout << "Cleanup took " << timer.elapsedMilliseconds() << "ms" << std::endl;

}


} /* namespace NetworKit */

