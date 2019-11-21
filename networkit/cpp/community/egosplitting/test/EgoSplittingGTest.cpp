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
#include "../../../io/EdgeListReader.h"
#include "../EgoSplitting.h"
#include "../../PLM.h"
#include "../../../io/METISGraphReader.h"

namespace NetworKit {

class EgoSplittingGTest : public testing::Test {
public:
	Graph testGraph;

	EgoSplittingGTest() {
//		EdgeListReader reader('\t', 0);
//		testGraph = reader.read("/home/armin/graphs/com-amazon.ungraph.txt");
//		EdgeListReader reader(' ', 0);
//		testGraph = reader.read("../input/lfr_om3.graph");
		METISGraphReader reader{};
//		testGraph = reader.read("../input/FB_Auburn71.graph");
		testGraph = reader.read("../input/FB_Caltech36.graph");
		testGraph.removeNode(0);
		testGraph.forNeighborsOf(1, [&](node v){
			testGraph.removeEdge(1, v);
		});
	}

	void testEgoSplitting(const std::map<std::string, std::string> &parameters) {
		for (std::size_t i = 0; i < 1; ++i) {
			//	std::function<Partition(const Graph &)> clusterAlgo = [](const Graph &G) {
			//		PLM plm(G, true, 1.0, "none");
			//		plm.run();
			//		return plm.getPartition();
			//	};
			//	EgoSplitting algo(G, clusterAlgo, clusterAlgo);
			bool parallelEgoNets = true;
			EgoSplitting algo(testGraph, parallelEgoNets);
			//	PLMFactory clusterFactory{true, 1.0, "none"};
			//	EgoSplitting algo(G, true, clusterFactory.getFunction(), clusterFactory.getFunction());
			algo.setParameters(parameters);
			algo.run();
			Cover cover = algo.getCover();

			std::cout << algo.timingsAsString() << std::endl;
			EXPECT_GE(cover.numberOfSubsets(), 5);
			for (auto size : cover.subsetSizes()) {
				EXPECT_GT(size, 4) << "discard communities with 4 or less nodes";
			}
		}
	}
};

TEST_F(EgoSplittingGTest, testEgoSplitting) {
//	Aux::Log::setLogLevel("INFO");
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

TEST_F(EgoSplittingGTest, testEgoSplittingCleanup) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Edges";
	parameters["Edges Score Strategy"] = "Edges pow 2 div Degree";
	parameters["Cleanup"] = "Yes";
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
