/*
 * LouvainMapEquation.h
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
#include "ClusteringFunctionFactory.h"

namespace NetworKit {

class LouvainMapEquation : public Algorithm {
public:
	explicit LouvainMapEquation(const Graph &graph, bool hierarchical = false, count maxIterations = 256);

	void run() override;

	/**
	 * Returns the result of the algorithm.
	 */
	const Partition& getPartition() const;

	std::string toString() const override;

private:
	const Graph &graph;
	bool hierarchical;
	count maxIterations;

	Partition partition;
	std::vector<double> clusterVolume;
	std::vector<double> clusterCut;
	double totalVolume;
	double totalCut;
	SparseVector<node> neighborClusterWeights;

	double fitnessChange(node, double degree, double loopWeight,
			     node currentCluster, node targetCluster,
			     double weightToTarget, double weightToCurrent);

        void moveNode(node u, double degree, double loopWeight, node currentCluster,
	              node targetCluster, double weightToTarget, double weightToCurrent);

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

class LouvainMapEquationFactory : public ClusteringFunctionFactory {
public:
	explicit LouvainMapEquationFactory(bool hierarchical = false, count maxIterations = 256);

	ClusteringFunction getFunction() const override;

private:
	bool hierarchical;
	count maxIterations;
};

}
