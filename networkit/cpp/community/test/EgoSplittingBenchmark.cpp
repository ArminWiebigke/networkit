/*
 * EgoSplittingBenchmark.cpp
 *
 * Created: 2019-10-18
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>

#include "../../auxiliary/Random.h"
#include "../PLP.h"
#include "../PLM.h"
#include "../Modularity.h"
#include "../EgoSplitting.h"
#include "../LPPotts.h"
#include "../../centrality/Betweenness.h"
#include "../../centrality/PageRank.h"
#include "../../auxiliary/Timer.h"
#include "../../structures/Partition.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"
#include "../../io/EdgeListReader.h"
#include "../../generators/ClusteredRandomGraphGenerator.h"

namespace NetworKit {

class EgoSplittingBenchmark: public testing::Test {

};

void benchEgoSplitting(const std::map<std::string, std::string>& parameters) {
//	EdgeListReader reader('\t', 0);
//	Graph G = reader.read("/home/armin/Code/graphs/com-amazon.ungraph.txt");
//	Cover C = CoverReader{}.read("/home/armin/Code/graphs/com-amazon.all.dedup.cmty.txt", G);
	EdgeListReader reader(' ', 0);
//	METISGraphReader reader{};
	Graph G = reader.read("../input/lfr_om3.graph");

	std::function<Partition(const Graph &)> clusterAlgo = [](const Graph &G) {
		PLM plm(G, true, 1.0, "none");
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

TEST_F(EgoSplittingBenchmark, benchSignificance) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Significance";
	parameters["Extend and Partition Iterations"] = "2";
	parameters["useSigMemo"] = "No";
	benchEgoSplitting(parameters);
}

TEST_F(EgoSplittingBenchmark, benchSignificanceMemoization) {
	std::map<std::string, std::string> parameters;
	parameters["Extend EgoNet Strategy"] = "Significance";
	parameters["Extend and Partition Iterations"] = "2";
	parameters["useSigMemo"] = "Yes";
	benchEgoSplitting(parameters);
}

} /* namespace NetworKit */
