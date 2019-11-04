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
#include "../structures/SparseVector.h"

namespace NetworKit {

class LocalMoveMapEquation : public Algorithm {
public:
	explicit LocalMoveMapEquation(Graph &graph, bool hierarchical = false, count maxIterations = 256);

	void run() override;

	/**
	 * Returns the result of the algorithm.
	 */
	Partition getPartition();

	std::string toString() const override;

private:
	Graph &graph;
	bool hierarchical;
	count maxIterations;

	Partition partition;
	std::vector<count> clusterVolume;
	std::vector<count> clusterCut;
	count totalVolume;
	count totalCut;
	SparseVector<node> neighborClusterWeights;

	double fitnessChange(node, count degree, count loopWeight, node currentCluster,
	                     node targetCluster, count weightToTarget,
	                     count weightToCurrent);

	void moveNode(node u, count degree, count loopWeight, node currentCluster,
	              node targetCluster, count weightToTarget, edgeid weightToCurrent);

#ifndef NDEBUG
	long double sumPLogPwAlpha = 0;
	long double sumPLogPClusterCut = 0;
	long double sumPLogPClusterCutPlusVol = 0;
	double plogpRel(count w);
	void updatePLogPSums();
	double mapEquation();
#endif

	bool tryLocalMove(node u);

	void calculateClusterCutAndVolume();

	void runHierarchical();
};

}

