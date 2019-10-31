/*
 * LocalMoveMapEquation.h
 *
 * Created on: 2019-01-28
 * Author: Armin Wiebigke
 *         Michael Hamann
 */

#include <vector>
#include <cstddef>
#include <algorithm>
#include <random>
#include <cmath>

#include "../Globals.h"
#include "../graph/Graph.h"
#include "../base/Algorithm.h"
#include "../structures/Partition.h"

namespace NetworKit {

class LocalMoveMapEquation : public Algorithm {
public:
	LocalMoveMapEquation(Graph &graph, bool hierarchical = false);

	LocalMoveMapEquation(Graph &graph, bool hierarchical,
	                     long double sum_p_log_p_w_alpha);

	void run() override;

	Partition getPartition();

	std::string toString() const override;

private:
	Graph &graph;
	bool hierarchical;
	long double sum_p_log_p_w_alpha;
	count max_rounds = 256;

	Partition partition;
	std::vector<count> clusterVolume;
	std::vector<count> clusterCut;
	count totalVolume;
	count totalCut;

	double fitnessChange(node, count degree, count loopWeight, node currentCluster,
	                     node targetCluster, count weightToTarget,
	                     count weightToCurrent);

	void moveNode(node u, count degree, count loopWeight, node currentCluster,
	              node targetCluster, count weightToTarget, edgeid weightToCurrent);

#ifndef NDEBUG
	long double sum_p_log_p_cluster_cut = 0;
	long double sum_p_log_p_cluster_cut_plus_vol = 0;
	double plogp_rel(count w);
	void update_p_log_p_sums();
	double map_equation();
#endif
};

}

