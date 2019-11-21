/*
 * LouvainMapEquation.cpp
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

#include "LouvainMapEquation.h"
#include "../auxiliary/SignalHandling.h"
#include "../graph/Graph.h"
#include "../structures/Partition.h"
#include "../coarsening/ParallelPartitionCoarsening.h"

namespace NetworKit {

LouvainMapEquation::LouvainMapEquation(const Graph &graph, bool hierarchical, count maxIterations)
		: graph(graph), hierarchical(hierarchical), maxIterations(maxIterations),
		  clusterVolume(graph.upperNodeIdBound()),
		  clusterCut(graph.upperNodeIdBound()), totalVolume(0), totalCut(0), partition(graph.upperNodeIdBound()),
		  neighborClusterWeights(graph.upperNodeIdBound(), none) {
}

void LouvainMapEquation::run() {
	if (hasRun)
		throw std::runtime_error("Algorithm was already run!");
	Aux::SignalHandler handler;
	partition.allToSingletons();
	for (node u = 0; u < graph.upperNodeIdBound(); ++u) {
		if (!graph.hasNode(u)) {
			partition.remove(u);
		}
	}
	handler.assureRunning();

	calculateClusterCutAndVolume();

#ifndef NDEBUG
	updatePLogPSums();

	if (sumPLogPwAlpha == 0) {
		sumPLogPwAlpha = std::accumulate<decltype(clusterCut.begin()), long double>(
				clusterVolume.begin(), clusterVolume.end(), .0, [&](long double sum, count vol) {
					return sum + plogpRel(vol);
				});
	}
#endif

	for (count iteration = 0; iteration < maxIterations; ++iteration) {
		handler.assureRunning();
		bool anyMoved = false;
		count nodesMoved = 0;
		INFO("\nIteration ", iteration);
#ifndef NDEBUG
		INFO("Map equation is ", mapEquation());
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
	handler.assureRunning();
	if (hierarchical && partition.numberOfSubsets() < graph.numberOfNodes()) {
		runHierarchical();
	}

	hasRun = true;
}

void LouvainMapEquation::calculateClusterCutAndVolume() {
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

bool LouvainMapEquation::tryLocalMove(node u) {
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

void LouvainMapEquation::runHierarchical() {
	assert(partition.numberOfSubsets() < partition.numberOfElements());
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

	LouvainMapEquation recursion(metaGraph, true, maxIterations);
	recursion.run();
	Partition metaPartition = recursion.getPartition();

	graph.forNodes([&](node u) {
		partition[u] = metaPartition[fineToCoarseMapping[u]];
	});
}

Partition LouvainMapEquation::getPartition() {
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
LouvainMapEquation::fitnessChange(node, count degree, count loopWeight, node currentCluster, node targetCluster,
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

void LouvainMapEquation::moveNode(node u, count degree, count loopWeight, node currentCluster,
                                  node targetCluster, count weightToTarget,
                                  count weightToCurrent) {
#ifndef NDEBUG
	long double oldVal = mapEquation();
	assert(oldVal > 0);
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
	updatePLogPSums();
	long double newVal = mapEquation();
	assert(newVal > 0);
//	std::cout << "Move node " << u << " from cluster " << currentCluster << " to " << targetCluster << std::endl;
//	std::cout << "Old: " << oldVal << ", fitnessDiff: " << fitnessDiff << " new: " << newVal << std::endl;
//	std::cout << "After update: " << mapEquation() << std::endl;
	assert(std::accumulate(clusterCut.begin(), clusterCut.end(), 0ull) == totalCut);
#endif
}

std::string LouvainMapEquation::toString() const {
	return "LouvainMapEquation";
}

#ifndef NDEBUG

double LouvainMapEquation::plogpRel(count w) {
	if (w > 0) {
		double p = static_cast<double>(w) / totalVolume;
		return p * log(p);
	}
	return 0;
}

void LouvainMapEquation::updatePLogPSums() {
	sumPLogPClusterCut = std::accumulate<decltype(clusterCut.begin()), long double>(
			clusterCut.begin(), clusterCut.end(), .0, [&](long double sum, count cut) {
				return sum + plogpRel(cut);
			});

	sumPLogPClusterCutPlusVol = 0;
	for (index i = 0; i < clusterCut.size(); ++i) {
		sumPLogPClusterCutPlusVol += plogpRel(clusterCut[i] + clusterVolume[i]);
	}
}

double LouvainMapEquation::mapEquation() {
	return plogpRel(totalCut) - 2 * sumPLogPClusterCut + sumPLogPClusterCutPlusVol -
	       sumPLogPwAlpha;
}

#endif

LouvainMapEquationFactory::LouvainMapEquationFactory(bool hierarchical, count maxIterations)
		: hierarchical(hierarchical), maxIterations(maxIterations) {
}

ClusteringFunction LouvainMapEquationFactory::getFunction() const {
	bool hiearchicalCopy = hierarchical;
	count maxIterationsCopy = maxIterations;
	return [hiearchicalCopy, maxIterationsCopy](const Graph &graph) {
		LouvainMapEquation algo(graph, hiearchicalCopy, maxIterationsCopy);
		algo.run();
		return algo.getPartition();
	};
}

}

