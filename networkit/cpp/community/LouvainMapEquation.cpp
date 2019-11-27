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
#include "../auxiliary/Timer.h"

namespace NetworKit {

LouvainMapEquation::LouvainMapEquation(const Graph &graph, bool hierarchical, count maxIterations, double additionalCut = 0.0, double additionalVolume = 0.0)
		: graph(graph), hierarchical(hierarchical), maxIterations(maxIterations),
		  clusterCut(graph.upperNodeIdBound()),
		  clusterVolume(graph.upperNodeIdBound()),
		  additionalCut(additionalCut),
		  additionalVolume(additionalVolume),
		  totalCut(additionalCut),
		  totalVolume(additionalVolume),
		  partition(graph.upperNodeIdBound()),
		  nextPartition(graph.upperNodeIdBound())
{

}

void LouvainMapEquation::run() {
	if (hasRun)
		throw std::runtime_error("Algorithm was already run!");
	
	Aux::SignalHandler handler;

	Aux::Timer timer;
	timer.start();
	
	partition.allToSingletons();
	nextPartition.allToSingletons();
	
	if (graph.numberOfNodes() != graph.upperNodeIdBound()) {
		for (node u = 0; u < graph.upperNodeIdBound(); ++u) {
			if (!graph.hasNode(u)) {
				partition.remove(u);
				nextPartition.remove(u);
			}
		}
	}
	handler.assureRunning();

	calculateClusterCutAndVolume();
	timer.stop();
	DEBUG("init ", timer.elapsedMilliseconds(), " ms");


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
	std::vector< SparseVector<double> > ets_neighborClusterWeights(Aux::getMaxNumberOfThreads());
	
	// could abuse neighborClusterWeights
	std::vector< std::vector<double> > ets_volumeUpdates(Aux::getMaxNumberOfThreads());
	std::vector< std::vector<double> > ets_cutUpdates(Aux::getMaxNumberOfThreads());
	
	for (count iteration = 0; iteration < maxIterations; ++iteration) {
		handler.assureRunning();

		DEBUG("Iteration ", iteration);
#ifndef NDEBUG
		DEBUG("Map equation is ", mapEquation());
#endif

		timer.start();
		std::shuffle(nodes.begin(), nodes.end(), Aux::Random::getURNG());
		timer.stop();
		DEBUG("shuffle ", timer.elapsedMilliseconds(), " ms");

		// chunks = fixed number of nodes, or degree sum ?
		const size_t chunkSize = 420;
		const size_t numberOfChunks = 1 + (graph.numberOfNodes() / chunkSize);
		std::vector<size_t> chunks;//chunks(numberOfChunks, 0);
		chunks.push_back(0);
		chunks.push_back(10);
		chunks.push_back(20);
		
		count numberOfNodesMoved = 0;
		
		timer.start();
		#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			
			SparseVector<double>& neighborClusterWeights = ets_neighborClusterWeights[tid];
			neighborClusterWeights.resize(graph.upperNodeIdBound(), 0.0);
			
			std::vector<double>& volumeUpdates = ets_volumeUpdates[tid];
			volumeUpdates.resize(graph.upperNodeIdBound(), 0.0);	// allocation happens once, not every iteration. and then it's local to the socket
			std::vector<double>& cutUpdates = ets_cutUpdates[tid];
			cutUpdates.resize(graph.upperNodeIdBound(), 0.0);
			
			std::vector<node> movedNodes;
			movedNodes.reserve(chunkSize);
			
			for (size_t i = 0; i < chunks.size() - 1; ++i) {
				const size_t firstInvalid = chunks[i + 1];
				
				// find moves
				#pragma omp for		// this works :) I wonder if distributing this loop across cores incurs some strange overhead
				for (size_t j = chunks[i]; j < firstInvalid; ++j) {
					const node u = nodes[j];
					if (tryLocalMove(u, neighborClusterWeights)) {
						movedNodes.push_back(u);
					}
				}
				
				// aggregate cut and volume updates
				// option a) every core iterates over its own performed moves, option b) iterate over the entire round again and check if partition changed
				// option a) has better locality, option b) better load balancing
				for (const node u : movedNodes) {
					const index cu = partition[u], ncu = nextPartition[u];
					assert(cu != ncu);
					double volU = 0.0;
					graph.forEdgesOf(u, [&](node , node v, edgeweight w) {
						if (cu != partition[v]) {
							cutUpdates[cu] -= w;
						}
						if (ncu != nextPartition[v]) {
							cutUpdates[ncu] += w;
						}
						volU += w;
						if (u == v) {
							volU += w;
						}
					});
					
					volumeUpdates[cu] -= volU;
					volumeUpdates[ncu] += volU;
				}
				
				auto applyUpdate = [&](std::vector<double>& global, std::vector<double>& local, index i) {
					if (local[i] != 0.0) {
						#pragma omp atomic
						global[i] += local[i];
						local[i] = 0.0;
					}
				};
				
				// apply aggregated updates to cut, volume and partition
				for (const node u : movedNodes) {
					const index cu = partition[u], ncu = nextPartition[u];
					applyUpdate(clusterVolume, volumeUpdates, cu);
					applyUpdate(clusterVolume, volumeUpdates, ncu);
					applyUpdate(clusterCut, cutUpdates, cu);
					applyUpdate(clusterCut, cutUpdates, ncu);
					partition[u] = nextPartition[u];
				}


				#ifndef NDEBUG
					#pragma omp single
					checkUpdatedCutsAndVolumesAgainstRecomputation();
				#endif



				#pragma omp atomic
				numberOfNodesMoved += movedNodes.size();
				movedNodes.clear();
			}
		}
		timer.stop();

		DEBUG("move iteration ", iteration, " took ", " ms");

		DEBUG("Moved ", numberOfNodesMoved, " nodes");
		clusteringChanged |= numberOfNodesMoved > 0;
		if (numberOfNodesMoved == 0) {
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

bool LouvainMapEquation::tryLocalMove(node u, SparseVector<double>& neighborClusterWeights) {
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

	if (neighborClusterWeights.size() > 0) {
		// Calculate best cluster
		node targetCluster = currentCluster;
		double bestChange = fitnessChange(u, degree, loop, currentCluster, currentCluster, weightToCurrent, weightToCurrent);
		for (index neighborCluster : neighborClusterWeights.insertedIndexes()) {
			const double neighborClusterWeight = neighborClusterWeights[neighborCluster];
			const double change = fitnessChange(u, degree, loop, currentCluster, neighborCluster, neighborClusterWeight, weightToCurrent);
			if (change < bestChange || (change == bestChange && neighborCluster < targetCluster)) {
				bestChange = change;
				targetCluster = neighborCluster;
			}
		}
		
		neighborClusterWeights.reset();
		
		// dont' apply move yet. just save it
		if (targetCluster != currentCluster) {
			nextPartition[u] = targetCluster;
			return true;
		}
	}
	return false;
}

void LouvainMapEquation::runHierarchical() {
	assert(partition.numberOfSubsets() < partition.numberOfElements());
	INFO("Run hierarchical with ", partition.numberOfSubsets(), " clusters (from ", graph.numberOfNodes(), " nodes)");
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
		targetCutPlusVolumeNew = cutTarget + cutDifferenceTarget + volTarget + degree;
		targetCutPlusVolumeCurrent = cutTarget + volTarget;
	} else {
		totalCutNew = totalCut;
		targetClusterCutNew = cutTarget;
		targetClusterCutCurrent = cutTarget + cutDifferenceCurrent;
		targetCutPlusVolumeNew = cutTarget + volTarget;
		targetCutPlusVolumeCurrent = cutTarget + cutDifferenceCurrent + volTarget - degree;
	}

	auto normalizeAndPLogP = [&](double& x) {
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

	return totalCutNew + ((targetCutPlusVolumeNew - targetCutPlusVolumeCurrent) - (2 * (targetClusterCutNew - targetClusterCutCurrent)));
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

void checkUpdatedCutsAndVolumesAgainstRecomputation() const {
	std::vector<double> cut, vol;
}

#endif

LouvainMapEquationFactory::LouvainMapEquationFactory(bool hierarchical, count maxIterations)
		: hierarchical(hierarchical), maxIterations(maxIterations) {
}

ClusteringFunction LouvainMapEquationFactory::getFunction() const {
	bool hierarchicalCopy = hierarchical;
	count maxIterationsCopy = maxIterations;
	return [hierarchicalCopy, maxIterationsCopy](const Graph &graph) {
		LouvainMapEquation algo(graph, hierarchicalCopy, maxIterationsCopy);
		algo.run();
		return algo.getPartition();
	};
}

}

