/*
 * OslomGTest.cpp
 *
 * Created: 2019-06-06
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>

#include "../../generators/ClusteredRandomGraphGenerator.h"
#include "../../community/EgoSplitting.h"
#include "../../structures/Cover.h"
#include "../../io/EdgeListReader.h"
#include "../../io/CoverReader.h"
#include "../OslomCleanUp.h"
#include "../Stochastics.h"
#include "../LogTable.h"

namespace NetworKit {

class OslomGTest : public testing::Test {
};

TEST_F(OslomGTest, testOslomCleanUp) {
	for (int i = 0; i < 5; ++i) {
		ClusteredRandomGraphGenerator gen(100, 4, 0.5, 0.02);
		Graph G = gen.generate();
		//	EdgeListReader reader('\t', 0);
		//	Graph G = reader.read("/home/armin/Code/graphs/com-amazon.ungraph.txt");
		//	Cover C = CoverReader{}.read("/home/armin/Code/graphs/com-amazon.all.dedup.cmty.txt",
		//								 G);
		//	EdgeListReader reader(' ', 0);
		//	Graph G = reader.read("/home/armin/Code/graphs/email-Eu-core.txt");


		EgoSplitting algo(G);
		algo.run();
		Cover cover = algo.getCover();

		std::vector<std::string> arguments{"-simple_cleanup",
		                                   "-merge_discarded", "-discard_max_extend_groups",
		                                   "-max_extend", "2",
		                                   "-cup_runs", "1",};
		OslomCleanUp cleanUp(G, cover, arguments);
		cleanUp.run();
		Cover cleanedCover = cleanUp.getCover();

		EXPECT_TRUE(cleanedCover.numberOfSubsets() <= cover.numberOfSubsets());
	}
}

TEST_F(OslomGTest, testOslomStochastics) {
	Stochastics::init();
	LogFactTable *log_table = LogFactTable::get_instance();
	log_table->set(100000);

	int avg_k = 50;
	int n = 1000;
	int external_stubs = avg_k * n - 3000;
	int group_out_stubs = 2000;
	int node_degree = avg_k;

	for (int k_in = 0; k_in < 15; ++k_in) {
		double r;
//		r = Stochastics::topological_05(k_in, group_out_stubs, external_stubs, avg_k);
		r = log_table->right_cumulative_function(node_degree, group_out_stubs, external_stubs,
		                                         k_in);
		double hyper = log_table->hyper(k_in, group_out_stubs, external_stubs, node_degree);
		double cum = Stochastics::order_statistics_left_cumulative(n - avg_k, n - avg_k, r);

		std::cout << k_in << ": r=" << r << ", hyper=" << hyper << ", cum=" << cum << std::endl;
//		std::cout << k_in << ": " << r << ", " << r * n << ", " << std::pow(1 - r, n) << std::endl;
	}
}

} /* namespace NetworKit */
