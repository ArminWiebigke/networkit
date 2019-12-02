/*
 * LouvainMapEquation.h
 *
 * Created on: 2019-01-28
 * Author: Armin Wiebigke
 *         Michael Hamann
 *         Lars Gottesbüren
 */

#include <vector>
#include <cstddef>
#include <algorithm>
#include <random>
#include <cmath>
#include <mutex>

#include "../Globals.h"
#include "../graph/Graph.h"
#include "../base/Algorithm.h"
#include "../structures/Partition.h"
#include "../structures/SparseVector.h"
#include "ClusteringFunctionFactory.h"

namespace NetworKit {


class LouvainMapEquation : public Algorithm {
public:
	explicit LouvainMapEquation(const Graph &graph, bool hierarchical = false, count maxIterations = 256, double additionalCut = 0.0, double additionalVolume = 0.0);

	void run() override;

	/**
	 * Returns the result of the algorithm.
	 */
	const Partition& getPartition() const;

	std::string toString() const override;

private:
	struct Move {
		node movedNode;
		double volume;
		index cacheID, originCluster, targetCluster;
		double cutUpdateToOriginCluster, cutUpdateToTargetCluster;
		
		Move(const node n = none, double vol = 0.0, index c = none, index cc = none, index tc = none, double cuptoc = 0.0, double cupttc = 0.0) :
				movedNode(n), volume(vol), cacheID(c), originCluster(cc), targetCluster(tc), cutUpdateToOriginCluster(cuptoc), cutUpdateToTargetCluster(cupttc) { }
	};
	
	struct NeighborInChunk {
		node neighbor;
		index oldCluster;
		double weightToNeighbor;
		NeighborInChunk(node n = none, index oc = none, double wtn = 0.0) : neighbor(n), oldCluster(oc), weightToNeighbor(wtn) { }
	};
	
	static_assert(std::is_trivially_destructible<Move>::value);
	static_assert(std::is_trivially_destructible<NeighborInChunk>::value);
	
	using NeighborCache = std::vector<NeighborInChunk>;
	using NeighborCaches = std::vector<NeighborCache>;
	
	
	const Graph& graph;
	bool hierarchical;
	count maxIterations;

	Partition partition, nextPartition;
	std::vector<double> clusterCut, clusterVolume;
	const double additionalCut, additionalVolume;		// only for debug mode to compare updated cuts/volumes against recomputed cuts/volumes
	double totalCut, totalVolume;
	
	
	double fitnessChange(node, double degree, double loopWeight,
			     node currentCluster, node targetCluster,
			     double weightToTarget, double weightToCurrent);

#ifndef NDEBUG
	long double sumPLogPwAlpha = 0;
	long double sumPLogPClusterCut = 0;
	long double sumPLogPClusterCutPlusVol = 0;
	double plogpRel(count w);
	void updatePLogPSums();
	double mapEquation();
	void checkUpdatedCutsAndVolumesAgainstRecomputation();
#endif

	bool tryLocalMove(node u, index& cacheID, std::vector<NeighborInChunk>& cachedNeighbors, SparseVector<double>& neighborClusterWeights, std::vector<Move>& moves, std::vector<bool>& isNodeInCurrentChunk);
	
	void aggregateAndApplyCutAndVolumeUpdates(std::vector<Move>& moves, NeighborCaches& neighborCaches);

	void calculateInitialClusterCutAndVolume();

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

