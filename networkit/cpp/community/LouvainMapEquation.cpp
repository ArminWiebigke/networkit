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
#include <atomic>
#include <omp.h>

#include "LouvainMapEquation.h"
#include "../auxiliary/SignalHandling.h"
#include "../graph/Graph.h"
#include "../structures/Partition.h"
#include "../coarsening/ParallelPartitionCoarsening.h"
#include "../auxiliary/Parallelism.h"

namespace NetworKit {

LouvainMapEquation::LouvainMapEquation(const Graph &graph, bool hierarchical, count maxIterations)
		: graph(graph), hierarchical(hierarchical), maxIterations(maxIterations),
		  clusterVolume(graph.upperNodeIdBound()),
		  clusterCut(graph.upperNodeIdBound()), totalVolume(0), totalCut(0), partition(graph.upperNodeIdBound()),
		  locks(graph.numberOfNodes()) {
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

	bool clusteringChanged = false;
	std::vector<node> nodes = graph.nodes();
	std::vector< SparseVector<node> > ets_neighborClusterWeights(Aux::getMaxNumberOfThreads(), SparseVector<node>(graph.upperNodeIdBound(), 0.0));
	
	for (count iteration = 0; iteration < maxIterations; ++iteration) {
		handler.assureRunning();
		
		INFO("Iteration ", iteration);
#ifndef NDEBUG
		INFO("Map equation is ", mapEquation());
#endif

		std::shuffle(nodes.begin(), nodes.end(), Aux::Random::getURNG());
		std::vector<count> ets_nodesMoved(Aux::getMaxNumberOfThreads(), 0);
		
		#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			count& nodesMoved = ets_nodesMoved[tid];
			SparseVector<node>& neighborClusterWeights = ets_neighborClusterWeights[tid];

			
			#pragma omp for
			for (size_t i = 0; i < nodes.size(); ++i) {
				if (tryLocalMove(nodes[i], neighborClusterWeights)) {
					nodesMoved += 1;
				}
			}
		}

		count nodesMoved = 0;
		for (count x : ets_nodesMoved) {
			nodesMoved += x;
		}
		
		INFO("Moved ", nodesMoved, " nodes");
		clusteringChanged |= nodesMoved > 0;
		if (nodesMoved == 0) {
			break;
		}
		
		// partition.compact(true);		// try compacting to increase locality
	}

	handler.assureRunning();
	if (hierarchical && clusteringChanged) {
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

bool LouvainMapEquation::tryLocalMove(node u, SparseVector<node>& neighborClusterWeights) {
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
		double bestChange = fitnessChange(u, degree, loop, currentCluster, currentCluster, weightToCurrent, weightToCurrent, totalCut);
		for (index neighborCluster : neighborClusterWeights.insertedIndexes()) {
			const double neighborClusterWeight = neighborClusterWeights[neighborCluster];
			const double change = fitnessChange(u, degree, loop, currentCluster, neighborCluster, neighborClusterWeight, weightToCurrent, totalCut);
			if (change < bestChange || (change == bestChange && neighborCluster < targetCluster)) {
				bestChange = change;
				targetCluster = neighborCluster;
			}
		}

		// Move node to best cluster
		if (targetCluster != currentCluster) {
			moved = moveNode(u, degree, loop, currentCluster, targetCluster, neighborClusterWeights[targetCluster], weightToCurrent);
		}

		neighborClusterWeights.reset();
	}
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
                                  double weightToTarget, double weightToCurrent,
								  const double totalCutCurrently /* copy of totalCut so it at least stays consistent for the fitnessChange calculations */) {
	
	const double cutTarget = clusterCut[targetCluster];
	const double volTarget = clusterVolume[targetCluster];
	const double cutDifferenceCurrent = 2 * weightToCurrent - degree + 2 * loopWeight;
	double totalCutNew, targetClusterCutNew, targetClusterCutCurrent, targetCutPlusVolumeNew, targetCutPlusVolumeCurrent;
	if (currentCluster != targetCluster) {
		double cutDifferenceTarget = degree - 2 * weightToTarget - 2 * loopWeight;

		totalCutNew = totalCutCurrently + cutDifferenceCurrent + cutDifferenceTarget;
		targetClusterCutNew = cutTarget + cutDifferenceTarget;
		targetClusterCutCurrent = cutTarget;
		targetCutPlusVolumeNew = cutTarget + cutDifferenceTarget +
                                         volTarget + degree;
		targetCutPlusVolumeCurrent = cutTarget + volTarget;
	} else {
		totalCutNew = totalCutCurrently;
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

bool LouvainMapEquation::moveNode(node u, double degree, double loopWeight, node currentCluster,
                                  node targetCluster, double weightToTarget,
                                  double weightToCurrent) {

	// lock currentCluster and targetCluster
	locks[ std::min(currentCluster, targetCluster) ].lock();
	locks[ std::max(currentCluster, targetCluster) ].lock();
	
	// recompute weightToCurrent and weightToTarget
	weightToCurrent = 0;
	weightToTarget = 0;
	graph.forEdgesOf(u, [&](node, node v, edgeweight weight) {
		if (u != v) {
			if (partition[v] == currentCluster) {
				weightToCurrent += weight;
			} else if (partition[v] == targetCluster) {
				weightToTarget += weight;
			}
		}
	});
	
	// perform move
	bool moved = false;
	const double totalCutCurrently = totalCut;
	const double fitnessCurrent = fitnessChange(u, degree, loopWeight, currentCluster, currentCluster, weightToCurrent, weightToCurrent, totalCutCurrently);
	const double fitnessTarget = fitnessChange(u, degree, loopWeight, currentCluster, targetCluster, weightToTarget, weightToCurrent, totalCutCurrently);
	if (fitnessTarget < fitnessCurrent) {
		double cutDifferenceCurrent = 2 * weightToCurrent - degree + 2 * loopWeight;
		double cutDifferenceTarget = degree - 2 * weightToTarget - 2 * loopWeight;
		clusterCut[currentCluster] += cutDifferenceCurrent;
		clusterCut[targetCluster] += cutDifferenceTarget;
		
		clusterVolume[currentCluster] -= degree;
		clusterVolume[targetCluster] += degree;

		#pragma omp atomic
		totalCut += cutDifferenceCurrent + cutDifferenceTarget;
		
		partition.moveToSubset(targetCluster, u);
		moved = true;
	}
	
	// unlock clusters again
	locks[ std::max(currentCluster, targetCluster) ].unlock();
	locks[ std::min(currentCluster, targetCluster) ].unlock();
	
	return moved;
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

