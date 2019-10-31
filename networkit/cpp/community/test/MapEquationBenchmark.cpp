/*
 * MapEquationBenchmark.cpp
 *
 * Created on: 2019-10-31
 * Author: Armin Wiebigke
  */

#include <gtest/gtest.h>

#include "../LocalMoveMapEquation.h"
#include "../../io/METISGraphReader.h"
#include "../../generators/ClusteredRandomGraphGenerator.h"
#include "../egosplitting/EgoSplitting.h"
#include "../PLM.h"

namespace NetworKit {

class MapEquationBenchmark : public testing::Test {
};

TEST_F(MapEquationBenchmark, benchLarge) {
	Aux::Random::setSeed(2342556, false);
	ClusteredRandomGraphGenerator generator(5000, 200, 0.5, 0.002);
	Graph G = generator.generate();
	Partition groundTruth = generator.getCommunities();
	Aux::Timer timer{};
	timer.start();

	LocalMoveMapEquation mapequation(G);
	mapequation.run();
	auto partition = mapequation.getPartition();

	timer.stop();
	std::cout << mapequation.toString() << " took " << timer.elapsedMilliseconds() << "ms" << std::endl;
	std::cout << partition.numberOfSubsets() << " clusters" << std::endl;
}

TEST_F(MapEquationBenchmark, benchLargeHierachical) {
	Aux::Random::setSeed(2342556, false);
	ClusteredRandomGraphGenerator generator(5000, 200, 0.5, 0.002);
	Graph G = generator.generate();
	Partition groundTruth = generator.getCommunities();
	Aux::Timer timer{};
	timer.start();

	LocalMoveMapEquation mapequation(G, true);
	mapequation.run();
	auto partition = mapequation.getPartition();

	timer.stop();
	std::cout << mapequation.toString() << " took " << timer.elapsedMilliseconds() << "ms" << std::endl;
	std::cout << partition.numberOfSubsets() << " clusters" << std::endl;
}

}
