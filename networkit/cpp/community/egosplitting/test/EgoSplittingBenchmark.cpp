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

namespace NetworKit {

class EgoSplittingBenchmark : public testing::Test {
public:
	Graph testGraph;

	EgoSplittingBenchmark() {
		Aux::Log::setLogLevel("INFO");
		EdgeListReader reader('\t', 0);
//		testGraph = reader.read("/home/armin/graphs/com-amazon.ungraph.txt");
		testGraph = reader.read("/home/armin/graphs/com-lj.ungraph.txt");
//		EdgeListReader reader(' ', 0);
//		testGraph = reader.read("input/lfr_om3.graph");
//		METISGraphReader reader{};
//		testGraph = reader.read("input/FB_Auburn71.graph");
	}

	void benchEgoSplitting(const std::map<std::string, std::string> &parameters) {
		Aux::Random::setSeed(3450441, false);
		bool egoNetsParallel = false;
		EgoSplitting algo(testGraph, egoNetsParallel);
		algo.setParameters(parameters);

		algo.run();
		Cover cover = algo.getCover();

		std::cout << algo.timingsAsString() << std::endl;
	}
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
	parameters["maxEgoNetsPartitioned"] = "50000";
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
