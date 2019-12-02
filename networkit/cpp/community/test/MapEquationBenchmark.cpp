/*
 * MapEquationBenchmark.cpp
 *
 * Created on: 2019-10-31
 * Author: Armin Wiebigke
  */

#include <gtest/gtest.h>

#include "../LouvainMapEquation.h"
#include "../../io/METISGraphReader.h"
#include "../../generators/ClusteredRandomGraphGenerator.h"
#include "../egosplitting/EgoSplitting.h"
#include "../PLM.h"
#include "../../io/SNAPGraphReader.h"

namespace NetworKit {

class MapEquationBenchmark : public testing::Test {
public:
	void SetUp() {
		Aux::Random::setSeed(435913, false);
	}
};

TEST_F(MapEquationBenchmark, benchLarge) {
	ClusteredRandomGraphGenerator generator(5000, 200, 0.5, 0.002);
	Graph G = generator.generate();
	Partition groundTruth = generator.getCommunities();
	Aux::Timer timer{};
	timer.start();

	LouvainMapEquation mapequation(G, false);
	mapequation.run();
	auto partition = mapequation.getPartition();

	timer.stop();
	std::cout << mapequation.toString() << " took " << timer.elapsedMilliseconds() << "ms" << std::endl;
	std::cout << partition.numberOfSubsets() << " clusters" << std::endl;
}

TEST_F(MapEquationBenchmark, benchLargeHierachical) {
//	ClusteredRandomGraphGenerator generator(5000, 200, 0.5, 0.002);
//	Graph G = generator.generate();
	std::cout << "[INPUT] graph file path (edge list tab 0, like SNAP) > " << std::endl;
	std::string graphPath;
	std::getline(std::cin, graphPath);

	Aux::Random::setSeed(420, true);

	SNAPGraphReader reader;
	Graph G = reader.read(graphPath);
//	Cover C = CoverReader{}.read("/home/armin/graphs/com-amazon.all.dedup.cmty.txt", G);
//	Partition groundTruth = generator.getCommunities();
	Aux::Timer timer{};
	timer.start();

	LouvainMapEquation mapequation(G, true);
	mapequation.run();
	auto partition = mapequation.getPartition();

	timer.stop();
	std::cout << mapequation.toString() << " took " << timer.elapsedMilliseconds() << "ms" << std::endl;
	std::cout << partition.numberOfSubsets() << " clusters" << std::endl;

}

}
