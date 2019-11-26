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
		  neighborClusterWeights(graph.upperNodeIdBound(), 0.0) {
}

void LouvainMapEquation::run() {
	if (hasRun)
		throw std::runtime_error("Algorithm was already run!");
	Aux::SignalHandler handler;
	partition.allToSingletons();
	if (graph.numberOfNodes() != graph.upperNodeIdBound()) {
		for (node u = 0; u < graph.upperNodeIdBound(); ++u) {
			if (!graph.hasNode(u)) {
				partition.remove(u);
			}
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
		DEBUG("Iteration ", iteration);
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

		DEBUG("Moved ", nodesMoved, " nodes");
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
	double degree = 0;
	double loop = 0;
	double weightToCurrent = 0;
	const node currentCluster = partition[u];
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

	bool moved = false;
	if (neighborClusterWeights.size() > 0) {
		// Calculate best cluster
		node targetCluster = currentCluster;
		double bestChange = fitnessChange(u, degree, loop, currentCluster, currentCluster, weightToCurrent, weightToCurrent);
		for (index neighborCluster : neighborClusterWeights.insertedIndexes()) {
			const double neighborClusterWeight = neighborClusterWeights[neighborCluster];
			const double change = fitnessChange(u, degree, loop, currentCluster, neighborCluster, neighborClusterWeight,
						    weightToCurrent);
			if (change < bestChange || (change == bestChange && neighborCluster < targetCluster)) {
				bestChange = change;
				targetCluster = neighborCluster;
			}
		}

		// Move node to best cluster
		if (targetCluster != currentCluster) {
			moveNode(u, degree, loop, currentCluster, targetCluster, neighborClusterWeights[targetCluster],
				weightToCurrent);
			moved = true;
		}

		neighborClusterWeights.reset();
	}

	return moved;
}

void LouvainMapEquation::runHierarchical() {
	assert(partition.numberOfSubsets() < partition.numberOfElements());
	INFO("Run hierarchical with ", partition.numberOfSubsets(), " clusters (from ", graph.numberOfNodes(), " nodes)");
	// free some memory
	clusterVolume.clear();
	clusterVolume.shrink_to_fit();
	clusterCut.clear();
	clusterCut.shrink_to_fit();
	neighborClusterWeights.clear();

	ParallelPartitionCoarsening coarsening(graph, partition, true, graph.numberOfNodes() > 1e6);
	coarsening.run();
	const Graph& metaGraph = coarsening.getCoarseGraph();
	const auto& fineToCoarseMapping = coarsening.getFineToCoarseNodeMapping();

	LouvainMapEquation recursion(metaGraph, true, maxIterations);
	recursion.run();
	const Partition& metaPartition = recursion.getPartition();

	graph.forNodes([&](node u) {
		partition[u] = metaPartition[fineToCoarseMapping[u]];
	});
}

const Partition& LouvainMapEquation::getPartition() const {
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
LouvainMapEquation::fitnessChange(node, double degree, double loopWeight, node currentCluster, node targetCluster,
                                  double weightToTarget, double weightToCurrent) {
	const double cutTarget = clusterCut[targetCluster];
	const double volTarget = clusterVolume[targetCluster];
	const double cutDifferenceCurrent = 2 * weightToCurrent - degree + 2 * loopWeight;
	double totalCutNew, targetClusterCutNew, targetClusterCutCurrent, targetCutPlusVolumeNew, targetCutPlusVolumeCurrent;
	if (currentCluster != targetCluster) {
		double cutDifferenceTarget = degree - 2 * weightToTarget - 2 * loopWeight;

		totalCutNew = totalCut + cutDifferenceCurrent + cutDifferenceTarget;
		targetClusterCutNew = cutTarget + cutDifferenceTarget;
		targetClusterCutCurrent = cutTarget;
		targetCutPlusVolumeNew = cutTarget + cutDifferenceTarget +
                                         volTarget + degree;
		targetCutPlusVolumeCurrent = cutTarget + volTarget;
	} else {
		totalCutNew = totalCut;
		targetClusterCutNew = cutTarget;
		targetClusterCutCurrent = cutTarget + cutDifferenceCurrent;
		targetCutPlusVolumeNew = cutTarget + volTarget;
		targetCutPlusVolumeCurrent =
		    cutTarget + cutDifferenceCurrent + volTarget - degree;
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

void LouvainMapEquation::moveNode(node u, double degree, double loopWeight, node currentCluster,
                                  node targetCluster, double weightToTarget,
                                  double weightToCurrent) {
#ifndef NDEBUG
	long double oldVal = mapEquation();
	assert(oldVal > 0);
	double moveFitness = fitnessChange(u, degree, loopWeight, currentCluster, targetCluster, weightToTarget,
	                                   weightToCurrent);
	double stayFitness = fitnessChange(u, degree, loopWeight, currentCluster, currentCluster, weightToCurrent,
	                                   weightToCurrent);
	double fitnessDiff = moveFitness - stayFitness;
#endif

	double cutDifferenceCurrent = 2 * weightToCurrent - degree + 2 * loopWeight;
	double cutDifferenceTarget = degree - 2 * weightToTarget - 2 * loopWeight;

	assert(clusterCut[currentCluster] >= cutDifferenceCurrent * -1);
	assert(clusterCut[targetCluster] >= cutDifferenceTarget * -1);
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

