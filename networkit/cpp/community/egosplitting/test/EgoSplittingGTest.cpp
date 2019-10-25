/*
 * EgoSplittingGTest.cpp
 *
 * Created: 2019-10-15
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>
#include <functional>

#include "../../../auxiliary/Random.h"
#include "../../../generators/ClusteredRandomGraphGenerator.h"
#include "../../../graph/Graph.h"
#include "../../../io/EdgeListReader.h"
#include "../EgoSplitting.h"
#include "../../PLM.h"
#include "../../../io/METISGraphReader.h"

namespace NetworKit {

class EgoSplittingGTest : public testing::Test {
};

void testEgoSplitting(const std::map<std::string, std::string>& parameters) {
	METISGraphReader reader{};
	Graph G = reader.read("../input/lfr_small.graph");

	std::function<Partition(const Graph &)> clusterAlgo = [](const Graph &G) {
		PLM plm(G, true, 1.0, "none");
		plm.run();
		return plm.getPartition();
	};

	EgoSplitting algo(G, clusterAlgo, clusterAlgo);
	algo.setParameters(parameters);
	algo.run();
	Cover cover = algo.getCover();

	EXPECT_GE(cover.numberOfSubsets(), 5);
	for (auto size : cover.subsetSizes()) {
		EXPECT_GT(size, 4) << "discard communities with 4 or less nodes";
	}
}

TEST_F(EgoSplittingGTest, testEgoSplitting) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "None";
	testEgoSplitting(parameters);
}

TEST_F(EgoSplittingGTest, testEgoSplittingEdges) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Edges";
	parameters["Edges Score Strategy"] = "Edges pow 2 div Degree";
	testEgoSplitting(parameters);
}

TEST_F(EgoSplittingGTest, testEgoSplittingEdgesSignificance) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Edges";
	parameters["Edges Score Strategy"] = "Significance";
	testEgoSplitting(parameters);
}

TEST_F(EgoSplittingGTest, testEgoSplittingSignificance) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Significance";
	parameters["Extend and Partition Iterations"] = "2";
	testEgoSplitting(parameters);
}

TEST_F(EgoSplittingGTest, DISABLED_testEgoSplittingSignificanceMemoization) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Significance";
	parameters["Extend and Partition Iterations"] = "2";
	parameters["useSigMemo"] = "Yes";
	testEgoSplitting(parameters);
}

} /* namespace NetworKit */
