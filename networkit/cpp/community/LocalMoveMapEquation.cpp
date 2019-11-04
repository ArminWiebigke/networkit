/*
 * LocalMoveMapEquation.cpp
 *
 * Created on: 2019-01-28
 * Author: Armin Wiebigke
 *         Michael Hamann
 */

#include <unordered_map>
#include <cassert>
#include <chrono>
#include <cmath>
#include <random>

#include "LocalMoveMapEquation.h"
#include "../graph/Graph.h"
#include "../structures/Partition.h"
#include "../coarsening/ParallelPartitionCoarsening.h"

namespace NetworKit {

LocalMoveMapEquation::LocalMoveMapEquation(Graph &graph, bool hierarchical, count maxIterations)
		: graph(graph), hierarchical(hierarchical), maxIterations(maxIterations),
		  clusterVolume(graph.upperNodeIdBound()),
		  clusterCut(graph.upperNodeIdBound()), totalVolume(0), totalCut(0), partition(graph.upperNodeIdBound()),
		  neighborClusterWeights(graph.upperNodeIdBound(), none){
}

void LocalMoveMapEquation::run() {
	if (hasRun)
		throw std::runtime_error("Algorithm was already run!");
	partition.allToSingletons();

	calculateClusterCutAndVolume();

#ifndef NDEBUG
	update_p_log_p_sums();

	if (sum_p_log_p_w_alpha == 0) {
		sum_p_log_p_w_alpha = std::accumulate<decltype(clusterCut.begin()), long double>(
				clusterVolume.begin(), clusterVolume.end(), .0, [&](long double sum, count vol) {
					return sum + plogp_rel(vol);
				});
	}
#endif

	for (size_t iteration = 0; iteration < maxIterations; ++iteration) {
		bool anyMoved = false;
		count nodesMoved = 0;
		INFO("\nIteration ", iteration);
#ifndef NDEBUG
		INFO("Map equation is ", map_equation());
#endif

		graph.forNodesInRandomOrder([&](node u) {
			bool moved = tryLocalMove(u);
			if (moved) {
				anyMoved = true;
				++nodesMoved;
			}
		});

		INFO("Moved ", nodesMoved, " nodes");
		if (!anyMoved) {
			break;
		}
	}

	partition.compact();
	if (hierarchical && partition.numberOfSubsets() < graph.numberOfNodes()) {
		runHierarchical();
	}

	hasRun = true;
}

void LocalMoveMapEquation::calculateClusterCutAndVolume() {
	graph.forNodes([&](node u) {
		graph.forEdgesOf(u, [&](node, node v, edgeweight weight) {
			if (u != v) {
				clusterCut[u] += weight;
				totalCut += weight;
			} else {
				weight *= 2; // loop weight counts twice
			}
			clusterVolume[u] += weight;
			totalVolume += weight;
		});
	});
}

bool LocalMoveMapEquation::tryLocalMove(node u) {
	// Find neighbor clusters
	count degree = 0;
	count loop = 0;
	count weightToCurrent = 0;
	node currentCluster = partition[u];
	graph.forEdgesOf(u, [&](node, node v, edgeweight weight) {
		degree += weight;
		if (u != v) {
			node neighborCluster = partition[v];
			if (neighborCluster == currentCluster) {
				weightToCurrent += weight;
			} else {
				if (!neighborClusterWeights.indexIsUsed(neighborCluster))
					neighborClusterWeights.insert(neighborCluster, 0);
				neighborClusterWeights[neighborCluster] += weight;
			}
		} else {
			loop += weight;
			degree += weight;
		}
	});
	assert(degree == graph.weightedDegree(u));

	// Calculate best cluster
	node targetCluster = currentCluster;
	double bestChange = fitnessChange(u, degree, loop, currentCluster, currentCluster, weightToCurrent,
	                                  weightToCurrent);
	for (index neighborCluster : neighborClusterWeights.insertedIndexes()) {
		count neighborClusterWeight = neighborClusterWeights[neighborCluster];
		double change = fitnessChange(u, degree, loop, currentCluster, neighborCluster, neighborClusterWeight,
		                              weightToCurrent);
		if (change < bestChange || (change == bestChange && neighborCluster < targetCluster)) {
			bestChange = change;
			targetCluster = neighborCluster;
		}
	}

	// Move node to best cluster
	bool moved = false;
	if (targetCluster != currentCluster) {
		moveNode(u, degree, loop, currentCluster, targetCluster, neighborClusterWeights[targetCluster],
		         weightToCurrent);
		moved = true;
	}

	neighborClusterWeights.reset();
	return moved;
}

void LocalMoveMapEquation::runHierarchical() {
	INFO("Run hierarchical with ", partition.numberOfSubsets(), " nodes (", graph.numberOfNodes());
	// free some memory
	clusterVolume.clear();
	clusterVolume.shrink_to_fit();
	clusterCut.clear();
	clusterCut.shrink_to_fit();
	neighborClusterWeights.clear();

	ParallelPartitionCoarsening coarsening(graph, partition);
	coarsening.run();
	Graph metaGraph = coarsening.getCoarseGraph();
	auto fineToCoarseMapping = coarsening.getFineToCoarseNodeMapping();

	LocalMoveMapEquation recursion(metaGraph, true, maxIterations);
	recursion.run();
	Partition metaPartition = recursion.getPartition();

	graph.forNodes([&](node u) {
		partition[u] = metaPartition[fineToCoarseMapping[u]];
	});
}

Partition LocalMoveMapEquation::getPartition() {
	return partition;
}

/**
 * Calculate the change in the map equation if the node is moved from its current cluster to the target cluster.
 * To simplify the calculation, we remove terms that are constant for all target clusters. As a result, "moving" the
 * node to its current cluster gives a value != 0, although the complete map equation would not change.
 * @param degree
 * @param loopWeight
 * @param currentCluster
 * @param targetCluster
 * @param weightToTarget
 * @param weightToCurrent
 * @return
 */
double
LocalMoveMapEquation::fitnessChange(node, count degree, count loopWeight, node currentCluster, node targetCluster,
                                    count weightToTarget, count weightToCurrent) {
	count cutDifferenceCurrent = 2 * weightToCurrent - degree + 2 * loopWeight;
	double totalCutNew, targetClusterCutCurrent, targetClusterCutNew, targetCutPlusVolumeNew, targetCutPlusVolumeCurrent;
	if (currentCluster != targetCluster) {
		count cutDifferenceTarget = degree - 2 * weightToTarget - 2 * loopWeight;

		totalCutNew = static_cast<double>(totalCut + cutDifferenceCurrent + cutDifferenceTarget);
		targetClusterCutNew = static_cast<double>(clusterCut[targetCluster] + cutDifferenceTarget);
		targetClusterCutCurrent = static_cast<double>(clusterCut[targetCluster]);
		targetCutPlusVolumeNew = static_cast<double>(clusterCut[targetCluster] + cutDifferenceTarget +
		                                             clusterVolume[targetCluster] + degree);
		targetCutPlusVolumeCurrent = static_cast<double>(clusterCut[targetCluster] + clusterVolume[targetCluster]);
	} else {
		totalCutNew = static_cast<double>(totalCut);
		targetClusterCutNew = static_cast<double>(clusterCut[currentCluster]);
		targetClusterCutCurrent = static_cast<double>(clusterCut[currentCluster] + cutDifferenceCurrent);
		targetCutPlusVolumeNew = static_cast<double>(clusterCut[currentCluster] + clusterVolume[currentCluster]);
		targetCutPlusVolumeCurrent = static_cast<double>(clusterCut[currentCluster] + cutDifferenceCurrent +
		                                                 clusterVolume[currentCluster] - degree);
	}

	auto normalizeAndPLogP = [&](double &x) {
		x /= totalVolume;
		if (x > .0) {
			x *= std::log(x);
		}
	};

	normalizeAndPLogP(totalCutNew);
	normalizeAndPLogP(targetClusterCutNew);
	normalizeAndPLogP(targetClusterCutCurrent);
	normalizeAndPLogP(targetCutPlusVolumeNew);
	normalizeAndPLogP(targetCutPlusVolumeCurrent);

	return totalCutNew + ((targetCutPlusVolumeNew - targetCutPlusVolumeCurrent)
	                      - (2 * (targetClusterCutNew - targetClusterCutCurrent)));
}

void LocalMoveMapEquation::moveNode(node u, count degree, count loopWeight, node currentCluster,
                                    node targetCluster, count weightToTarget,
                                    count weightToCurrent) {
#ifndef NDEBUG
	long double old_val = map_equation();
	assert(old_val > 0);
	double moveFitness = fitnessChange(u, degree, loopWeight, currentCluster, targetCluster, weightToTarget,
									   weightToCurrent);
	double stayFitness = fitnessChange(u, degree, loopWeight, currentCluster, currentCluster, weightToCurrent,
									   weightToCurrent);
	double fitnessDiff = moveFitness - stayFitness;
#endif

	int64_t cutDifferenceCurrent = 2 * weightToCurrent - degree + 2 * loopWeight;
	int64_t cutDifferenceTarget = degree - 2 * weightToTarget - 2 * loopWeight;

	assert(static_cast<int64_t>(clusterCut[currentCluster]) >= cutDifferenceCurrent * -1);
	assert(static_cast<int64_t>(clusterCut[targetCluster]) >= cutDifferenceTarget * -1);
	clusterCut[currentCluster] += cutDifferenceCurrent;
	clusterCut[targetCluster] += cutDifferenceTarget;
	totalCut += cutDifferenceCurrent + cutDifferenceTarget;

	clusterVolume[currentCluster] -= degree;
	clusterVolume[targetCluster] += degree;

	partition.moveToSubset(targetCluster, u);

#ifndef NDEBUG
	update_p_log_p_sums();
	long double new_val = map_equation();
	assert(new_val > 0);
//	std::cout << "Move node " << u << " from cluster " << currentCluster << " to " << targetCluster << std::endl;
//	std::cout << "Old: " << old_val << ", fitnessDiff: " << fitnessDiff << " new: " << new_val << std::endl;
//	std::cout << "After update: " << map_equation() << std::endl;
	assert(std::accumulate(clusterCut.begin(), clusterCut.end(), 0ull) == totalCut);
#endif
}

std::string LocalMoveMapEquation::toString() const {
	return "LocalMoveMapEquation";
}

#ifndef NDEBUG

double LocalMoveMapEquation::plogp_rel(count w) {
	if (w > 0) {
		double p = static_cast<double>(w) / totalVolume;
		return p * log(p);
	}
	return 0;
}

void LocalMoveMapEquation::update_p_log_p_sums() {
	sum_p_log_p_cluster_cut = std::accumulate<decltype(clusterCut.begin()), long double>(
			clusterCut.begin(), clusterCut.end(), .0, [&](long double sum, count cut) {
				return sum + plogp_rel(cut);
			});

	sum_p_log_p_cluster_cut_plus_vol = 0;
	for (size_t i = 0; i < clusterCut.size(); ++i) {
		sum_p_log_p_cluster_cut_plus_vol += plogp_rel(clusterCut[i] + clusterVolume[i]);
	}
}

double LocalMoveMapEquation::map_equation() {
	return plogp_rel(totalCut) - 2 * sum_p_log_p_cluster_cut + sum_p_log_p_cluster_cut_plus_vol -
		   sum_p_log_p_w_alpha;
}

#endif

}
