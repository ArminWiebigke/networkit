/*
 * EgoSplittingGTest.cpp
 *
 * Created: 2019-10-15
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>
#include <functional>

#include "../../auxiliary/Random.h"
#include "../../generators/ClusteredRandomGraphGenerator.h"
#include "../../graph/Graph.h"
#include "../../io/EdgeListReader.h"
#include "../EgoSplitting.h"
#include "../PLM.h"

namespace NetworKit {

class EgoSplittingGTest : public testing::Test {
};

void testEgoSplitting(std::map<std::string, std::string> parameters) {
	Aux::Random::setSeed(234769, false);

//	EdgeListReader reader('\t', 0);
//	Graph G = reader.read("/home/armin/Code/graphs/com-amazon.ungraph.txt");
//	Cover C = CoverReader{}.read("/home/armin/Code/graphs/com-amazon.all.dedup.cmty.txt",
//								 G);
//	EdgeListReader reader(' ', 0);
//	Graph G = reader.read("../input/lfr_om3.graph");
	// TODO: small graph with overlapping communities
	ClusteredRandomGraphGenerator gen(100, 10, 0.5, 0.03);
	Graph G = gen.generate();
//	G.removeSelfLoops();
//	G.indexEdges();

	std::function<Partition(const Graph &)> clusterAlgo = [](const Graph &G) {
		PLM plm(G, false, 1.0, "none");
		plm.run();
		return plm.getPartition();
	};

	EgoSplitting algo(G, clusterAlgo, clusterAlgo);
	algo.setParameters(parameters);
	algo.run();
	Cover cover = algo.getCover();

	for (auto size : cover.subsetSizes()) {
		EXPECT_GT(size, 4) << "discard communities with 4 or less nodes";
	}
}

TEST_F(EgoSplittingGTest, testEgoSplittingEdges) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Edges";
	parameters["Edges Score Strategy"] = "Edges pow 2 div Degree";
	testEgoSplitting(parameters);
}

TEST_F(EgoSplittingGTest, testEgoSplittingSignificance) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Significance";
	parameters["Extend and Partition Iterations"] = "2";
	testEgoSplitting(parameters);
}

TEST_F(EgoSplittingGTest, testEgoSplittingSignificanceMemoization) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Significance";
	parameters["Extend and Partition Iterations"] = "2";
	parameters["useSigMemo"] = "Yes";
	testEgoSplitting(parameters);
}

} /* namespace NetworKit */
