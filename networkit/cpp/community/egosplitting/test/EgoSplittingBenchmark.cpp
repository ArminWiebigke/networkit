/*
 * EgoSplittingBenchmark.cpp
 *
 * Created: 2019-10-18
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>

#include "../../../auxiliary/Random.h"
#include "../../PLM.h"
#include "../../Modularity.h"
#include "../EgoSplitting.h"
#include "../../../structures/Partition.h"
#include "../../../io/METISGraphReader.h"
#include "../../../io/EdgeListReader.h"
#include "../../../io/SNAPGraphReader.h"

namespace NetworKit {

class EgoSplittingBenchmark : public testing::Test {
public:
	void SetUp() {
		Aux::Random::setSeed(435913, false);
	}

	EgoSplittingBenchmark() {
//		Aux::Log::setLogLevel("INFO");
//		EdgeListReader reader('\t', 0);
//		testGraph = reader.read("/home/armin/graphs/com-amazon.ungraph.txt");
//		EdgeListReader reader(' ', 0);
//		testGraph = reader.read("input/lfr_om3.graph");
		std::string graphPath;
		std::cout << "[INPUT] graph file path (edge list tab 0, like SNAP) > " << std::endl;
		std::getline(std::cin, graphPath);
		//EdgeListReader reader('\t', 0);
		SNAPGraphReader reader;
		testGraph = reader.read(graphPath);
	}

	void benchEgoSplitting(const std::map<std::string, std::string> &parameters) {
		Aux::Random::setSeed(3450441, false);
		bool egoNetsParallel = true;
		EgoSplitting algo(testGraph, egoNetsParallel);
		algo.setParameters(parameters);

		algo.run();
		Cover cover = algo.getCover();

		std::cout << algo.timingsAsString() << std::endl;
	}

	Graph testGraph;
};

TEST_F(EgoSplittingBenchmark, benchNoExtend) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "None";
	benchEgoSplitting(parameters);
}

TEST_F(EgoSplittingBenchmark, benchEdges) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Edges";
	parameters["Edges Score Strategy"] = "Edges pow 2 div Degree";
	benchEgoSplitting(parameters);
}

TEST_F(EgoSplittingBenchmark, benchEdgesNoMerge) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Edges";
	parameters["Edges Score Strategy"] = "Edges pow 2 div Degree";
	parameters["CleanupMerge"] = "No";
	benchEgoSplitting(parameters);
}

TEST_F(EgoSplittingBenchmark, benchSignificance) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Significance";
	parameters["Extend and Partition Iterations"] = "2";
	parameters["useSigMemo"] = "No";
	benchEgoSplitting(parameters);
}

TEST_F(EgoSplittingBenchmark, DISABLED_benchSignificanceMemoization) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Significance";
	parameters["Extend and Partition Iterations"] = "2";
	parameters["useSigMemo"] = "Yes";
	benchEgoSplitting(parameters);
}

} /* namespace NetworKit */
