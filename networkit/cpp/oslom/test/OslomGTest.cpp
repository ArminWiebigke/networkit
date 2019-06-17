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
	for (int i = 0; i < 1; ++i) {
		ClusteredRandomGraphGenerator gen(400, 10, 0.7, 0.05);
		Graph G = gen.generate();
		EdgeListReader reader('\t', 0);
//		Graph G = reader.read("/home/armin/Code/graphs/com-amazon.ungraph.txt");
//		Cover C = CoverReader{}.read("/home/armin/Code/graphs/com-amazon.all.dedup.cmty.txt",
//		                             G);
//		EdgeListReader reader(' ', 0);
//		Graph G = reader.read("/home/armin/Code/graphs/email-Eu-core.txt");
//		G.removeSelfLoops();

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
		count notEmptySubsets = 0;
		for (count s : cleanedCover.subsetSizes()) {
			notEmptySubsets += (s > 1);
		}
		EXPECT_GT(notEmptySubsets, 5);
	}
}

double calc_score(int node_degree, int k_in, int gr_out, int ext_stubs, int ext_nodes) {
	double boot_interval;
	double r_score = Stochastics::compute_simple_fitness(k_in, gr_out, ext_stubs,
	                                                     node_degree);
	double ordered_stat = Stochastics::order_statistics_left_cumulative(ext_nodes,
	                                                                    ext_nodes, r_score);
	return ordered_stat;
}

TEST_F(OslomGTest, testOslomStochastics) {
	int avg_k = 50;
	int n = 1000;
	int total_edges = avg_k * n;
	int node_degree = avg_k;
	int group_size = 4;
	int group_out_stubs = 40 * group_size;
	int group_internal_stubs = group_out_stubs / 2;
	int external_stubs = total_edges * 2 - group_out_stubs + group_internal_stubs;
	int external_nodes = n - group_size;
//	external_nodes = 10;

	Stochastics::init(total_edges);

	for (int k_in = 0; k_in < 15; ++k_in) {
//		double r;
//		r = Stochastics::topological_05(k_in, group_out_stubs, external_stubs, avg_k);
//		r = log_table->right_cumulative_function(node_degree, group_out_stubs, external_stubs,		                                         k_in);
//		double hyper = log_table->hypergeom_dist(k_in, group_out_stubs, external_stubs, node_degree);
//		double cum = Stochastics::order_statistics_left_cumulative(n - avg_k, n - avg_k, r);
//		std::cout << k_in << ": r=" << r << ", hyper=" << hypergeom_dist << ", cum=" << cum << std::endl;
////		std::cout << k_in << ": " << r << ", " << r * n << ", " << std::pow(1 - r, n) << std::endl;
		if (k_in > group_size)
			break;
		double boot_interval;
		double r_score = Stochastics::compute_simple_fitness(k_in, group_out_stubs, external_stubs,
		                                                     node_degree);
		double ordered_stat = Stochastics::order_statistics_left_cumulative(external_nodes,
		                                                                    external_nodes,
		                                                                    r_score);
//		if (ordered_stat < 0.1) {
			std::cout << k_in << ": fit=" << r_score << ", boot_interval=" << boot_interval
			          << ", ord_stat=" << ordered_stat << std::endl;
//			break;
//		}
//		std::cout << k_in << ": "
//		          << calc_score(node_degree, k_in, group_out_stubs, external_stubs, external_nodes)
//		          << std::endl;

	}
}


} /* namespace NetworKit */
