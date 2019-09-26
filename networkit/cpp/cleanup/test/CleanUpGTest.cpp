/*
 * OslomGTest.cpp
 *
 * Created: 2019-06-06
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>
#include <cmath>

#include "../../generators/ClusteredRandomGraphGenerator.h"
#include "../../community/EgoSplitting.h"
#include "../../structures/Cover.h"
#include "../../io/EdgeListReader.h"
#include "../../io/CoverReader.h"
#include "../SignificanceCommunityCleanUp.h"
#include "../StochasticDistribution.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {

class CleanupGTest : public testing::Test {
};

TEST_F(CleanupGTest, testCleanUp) {
	for (int i = 0; i < 3; ++i) {
		ClusteredRandomGraphGenerator gen(100, 10, 0.7, 0.05);
		Graph G = gen.generate();

//		EdgeListReader reader('\t', 0);
//		Graph G = reader.read("/home/armin/graphs/com-amazon.ungraph.txt");
//		Cover C = CoverReader{}.read("/home/armin/graphs/com-amazon.all.dedup.cmty.txt",
//		                             G);
//		EdgeListReader reader(' ', 0);
//		Graph G = reader.read("/home/armin/graphs/lfr_om3.txt");
//		Graph G = reader.read("/home/armin/graphs/email-Eu-core.txt");
//		G.removeSelfLoops();

//		METISGraphReader reader;
//		Graph G = reader.read("../input/jazz.graph");

		node isolatedNode = G.addNode();

		EgoSplitting algo(G);
		algo.run();
		Cover cover = algo.getCover();
		cover.addSubset({1});
		cover.addSubset({2, isolatedNode});

		SignificanceCommunityCleanUp cleanUp(G, cover, 0.1, 0.1, 0.5);
		cleanUp.run();
		Cover cleanedCover = cleanUp.getCover();

		std::cout << "Cleaned communities: " << cleanedCover.numberOfSubsets() << std::endl;
		EXPECT_TRUE(cleanedCover.numberOfSubsets() <= cover.numberOfSubsets());
		count notEmptyComms = 0;
		const std::vector<count> &comms = cleanedCover.subsetSizes();
		for (count s : comms) {
			notEmptyComms += (s > 1);
		}
		EXPECT_GE(notEmptyComms, 9);
		// Communities of size 1 should be discarded
		for (count s : comms) {
			if (s > 0)
				EXPECT_GT(s, 1);
		}
		std::set<index> badComm = {2, isolatedNode};
		for (const auto &comm : cleanedCover.getSubsets()) {
			EXPECT_NE(comm, badComm);
		}
	}
}


TEST_F(CleanupGTest, testBinomialCoeff) {
	StochasticDistribution stoch(10);

	EXPECT_NEAR(1, stoch.binomCoeff(7, 0), 1e-6);
	EXPECT_NEAR(7, stoch.binomCoeff(7, 1), 1e-6);
	EXPECT_NEAR(21, stoch.binomCoeff(7, 2), 1e-6);
	EXPECT_NEAR(35, stoch.binomCoeff(7, 3), 1e-6);
	EXPECT_NEAR(35, stoch.binomCoeff(7, 4), 1e-6);
	EXPECT_NEAR(21, stoch.binomCoeff(7, 5), 1e-6);
	EXPECT_NEAR(7, stoch.binomCoeff(7, 6), 1e-6);
	EXPECT_NEAR(1, stoch.binomCoeff(7, 7), 1e-6);
}

TEST_F(CleanupGTest, testBinomialDist) {
	StochasticDistribution stoch(10);

	EXPECT_NEAR(stoch.binomialDist(0.3, 7, 0), 0.0823543, 1e-7);
	EXPECT_NEAR(stoch.binomialDist(0.3, 7, 1), 0.247063, 1e-6);
	EXPECT_NEAR(stoch.binomialDist(0.3, 7, 2), 0.317652, 1e-6);
	EXPECT_NEAR(stoch.binomialDist(0.3, 7, 3), 0.226894, 1e-6);
	EXPECT_NEAR(stoch.binomialDist(0.3, 7, 4), 0.0972405, 1e-7);
	EXPECT_NEAR(stoch.binomialDist(0.3, 7, 5), 0.0250047, 1e-7);
	EXPECT_NEAR(stoch.binomialDist(0.3, 7, 6), 0.0035721, 1e-8);
	EXPECT_NEAR(stoch.binomialDist(0.3, 7, 7), 0.0002187, 1e-9);
}

TEST_F(CleanupGTest, testRightCumBinom) {
	StochasticDistribution stoch(10);

	EXPECT_NEAR(stoch.rightCumulativeBinomial(0.3, 5, 0), 1, 1e-6);
	EXPECT_NEAR(stoch.rightCumulativeBinomial(0.3, 5, 1), 0.83193, 1e-6);
	EXPECT_NEAR(stoch.rightCumulativeBinomial(0.3, 5, 2), 0.47178, 1e-6);
	EXPECT_NEAR(stoch.rightCumulativeBinomial(0.3, 5, 3), 0.16308, 1e-6);
	EXPECT_NEAR(stoch.rightCumulativeBinomial(0.3, 5, 4), 0.03078, 1e-7);
	EXPECT_NEAR(stoch.rightCumulativeBinomial(0.3, 5, 5), 0.00243, 1e-8);
}

TEST_F(CleanupGTest, testLeftCumBinom) {
	StochasticDistribution stoch(10);

	EXPECT_NEAR(stoch.leftCumulativeBinomial(0.3, 5, 0), 0.16807, 1e-6);
	EXPECT_NEAR(stoch.leftCumulativeBinomial(0.3, 5, 1), 0.52822, 1e-6);
	EXPECT_NEAR(stoch.leftCumulativeBinomial(0.3, 5, 2), 0.83692, 1e-6);
	EXPECT_NEAR(stoch.leftCumulativeBinomial(0.3, 5, 3), 0.96922, 1e-6);
	EXPECT_NEAR(stoch.leftCumulativeBinomial(0.3, 5, 4), 0.99757, 1e-6);
	EXPECT_NEAR(stoch.leftCumulativeBinomial(0.3, 5, 5), 1, 1e-6);
}

TEST_F(CleanupGTest, testHyperDist) {
	StochasticDistribution stoch(10);

	EXPECT_NEAR(stoch.hypergeometricDist(10, 5, 5, 0), 0.003968254, 1e-8);
	EXPECT_NEAR(stoch.hypergeometricDist(10, 5, 5, 1), 0.099206349, 1e-7);
	EXPECT_NEAR(stoch.hypergeometricDist(10, 5, 5, 2), 0.396825397, 1e-6);
	EXPECT_NEAR(stoch.hypergeometricDist(10, 5, 5, 3), 0.396825397, 1e-6);
	EXPECT_NEAR(stoch.hypergeometricDist(10, 5, 5, 4), 0.099206349, 1e-7);
	EXPECT_NEAR(stoch.hypergeometricDist(10, 5, 5, 5), 0.003968254, 1e-8);
}

TEST_F(CleanupGTest, testRightCumHyper) {
	StochasticDistribution stoch(10);

	EXPECT_NEAR(stoch.rightCumulativeHyper(10, 5, 5, 0), 1, 1e-6);
	EXPECT_NEAR(stoch.rightCumulativeHyper(10, 5, 5, 1), 0.996031746, 1e-6);
	EXPECT_NEAR(stoch.rightCumulativeHyper(10, 5, 5, 2), 0.896825397, 1e-6);
	EXPECT_NEAR(stoch.rightCumulativeHyper(10, 5, 5, 3), 0.5, 1e-6);
	EXPECT_NEAR(stoch.rightCumulativeHyper(10, 5, 5, 4), 0.103174603, 1e-6);
	EXPECT_NEAR(stoch.rightCumulativeHyper(10, 5, 5, 5), 0.003968254, 1e-8);
}

TEST_F(CleanupGTest, testLeftCumHyper) {
	StochasticDistribution stoch(10);

	EXPECT_NEAR(stoch.leftCumulativeHyper(10, 5, 5, 0), 0.003968254, 1e-8);
	EXPECT_NEAR(stoch.leftCumulativeHyper(10, 5, 5, 1), 0.103174603, 1e-6);
	EXPECT_NEAR(stoch.leftCumulativeHyper(10, 5, 5, 2), 0.5, 1e-6);
	EXPECT_NEAR(stoch.leftCumulativeHyper(10, 5, 5, 3), 0.896825397, 1e-6);
	EXPECT_NEAR(stoch.leftCumulativeHyper(10, 5, 5, 4), 0.996031746, 1e-6);
	EXPECT_NEAR(stoch.leftCumulativeHyper(10, 5, 5, 5), 1, 1e-6);
}

TEST_F(CleanupGTest, testStochasticDist) {
	StochasticDistribution stoch(10);
	count kTotal = 10;
	count kIn = 3;
	count cOut = 20;
	count extStubs = 100;
	count M = extStubs - kTotal;
	auto fct = [](count x) {
		double f = 1.0;
		for (count i = 2; i <= x; ++i)
			f *= i;
		return f;
	};
	auto p = [&](count kIn) {
		count kOut = kTotal - kIn;
		count MInEdges = (M - cOut - kOut + kIn) / 2;
		double top = std::pow(2, - (double) kIn);
		double divisor = fct(kOut) * fct(kIn) * fct(cOut - kIn) * fct(MInEdges);
		return top / divisor;
	};
	double probabilitySum = 0.0;
	double cumulativeProbSum = 0.0;
	for (count x = 0; x < kTotal; ++x) {
		double probability = p(x);
		probabilitySum += probability;
		if (x >= kIn)
			cumulativeProbSum += probability;
	}
	double exactProbCorrect = p(kIn) / probabilitySum;
	double cumulativeProbCorrect = cumulativeProbSum / probabilitySum;

	double exactProb, cumulativeProb;
	std::tie(exactProb, cumulativeProb) = stoch.rightCumulativeStochastic(
			kTotal, kIn, cOut, extStubs);

	EXPECT_NEAR(exactProb, exactProbCorrect, 1e-6);
	EXPECT_NEAR(exactProb / exactProbCorrect, 1.0, 1e-6);
	EXPECT_NEAR(cumulativeProb, cumulativeProbCorrect, 1e-6);
	EXPECT_NEAR(cumulativeProb / cumulativeProbCorrect, 1.0,1e-6);

}

} /* namespace NetworKit */
