/*
 * LouvainMapEquation.cpp
 *
 * Created on: 2019-01-28
 * Author: Armin Wiebigke
 *         Michael Hamann
 *         Lars Gottesbüren
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
#include "../auxiliary/Parallel.h"
#include "../auxiliary/Timer.h"

namespace NetworKit {

LouvainMapEquation::LouvainMapEquation(const Graph &graph, bool hierarchical, count maxIterations, double additionalCut, double additionalVolume)
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

	calculateInitialClusterCutAndVolume();
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
	
	std::vector< std::vector<double> > ets_volumeUpdates(Aux::getMaxNumberOfThreads());
	std::vector< SparseVector<double> > ets_cutUpdates(Aux::getMaxNumberOfThreads(), SparseVector<double>(0, std::numeric_limits<double>::infinity()));
	
	
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

		// chunkBorders = fixed number of nodes, or degree sum ?
		const size_t numberOfChunks = 5;
		const size_t chunkSize = std::max(size_t(1), nodes.size() / numberOfChunks);
		std::vector<size_t> chunkBorders = Aux::Parallel::Chunking::getChunkBorders(nodes.size(), numberOfChunks);
		
		count numberOfNodesMoved = 0;
		
		timer.start();
		#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			
			SparseVector<double>& neighborClusterWeights = ets_neighborClusterWeights[tid];
			neighborClusterWeights.resize(graph.upperNodeIdBound(), 0.0);
			std::vector<double>& volumeUpdates = ets_volumeUpdates[tid];
			volumeUpdates.resize(graph.upperNodeIdBound(), 0.0);	// allocation happens once. not every iteration. and then it's local to the socket. with thread pinning it remains local.
			SparseVector<double>& cutUpdates = ets_cutUpdates[tid];
			cutUpdates.resize(graph.upperNodeIdBound(), std::numeric_limits<double>::infinity());
			
			std::vector<node> movedNodes;
			movedNodes.reserve(chunkSize);

			for (size_t i = 0; i < chunkBorders.size() - 1; ++i) {
				
				// find moves
				const size_t firstInvalid = chunkBorders[i + 1];
				#pragma omp for
				for (size_t j = chunkBorders[i]; j < firstInvalid; ++j) {
					const node u = nodes[j];
					if (tryLocalMove(u, neighborClusterWeights)) {
						movedNodes.push_back(u);
					}
				}
				// implicit barrier at the end of the for-loop
				//#pragma omp barrier
				
				aggregateAndApplyCutAndVolumeUpdates(movedNodes, cutUpdates, volumeUpdates);
				#pragma omp barrier

				for (const node u : movedNodes) {
					assert(partition[u] != nextPartition[u]);
					partition[u] = nextPartition[u];
				}
				#pragma omp barrier
				
				
				#ifndef NDEBUG
					#pragma omp single
					{
						checkUpdatedCutsAndVolumesAgainstRecomputation();
					}
				#endif
				

				#pragma omp atomic update
				numberOfNodesMoved += movedNodes.size();
				movedNodes.clear();
				
			}
		}
		timer.stop();

		clusteringChanged |= numberOfNodesMoved > 0;
		if (numberOfNodesMoved == 0) {
			break;
		}
		
		DEBUG("Move iteration ", iteration, " took ", timer.elapsedMilliseconds(), " ms. Moved ", numberOfNodesMoved, " nodes");
	}

	handler.assureRunning();
	if (hierarchical && clusteringChanged) {
		runHierarchical();
	}
	hasRun = true;
}

void LouvainMapEquation::aggregateAndApplyCutAndVolumeUpdates(std::vector<node>& movedNodes, SparseVector<double>& cutUpdates, std::vector<double>& volumeUpdates) {
	// aggregate cut and volume updates
	// option a) every core iterates over its own performed moves, option b) iterate over the entire round in parallel again and check if partition changed
	// option a) has better locality, option b) has better load balancing but incurs an additional synchronization
	
	assert(cutUpdates.isClean());
	assert(std::all_of(volumeUpdates.begin(), volumeUpdates.end(), [](const double x) { return x == 0.0; }));
	
	double totalCutUpdate = 0.0;
	for (const node u : movedNodes) {
		const index cu = partition[u], ncu = nextPartition[u];
		assert(cu != ncu);
		double cutUpdateCU = 0.0, cutUpdateNCU = 0.0, otherCuts = 0.0, volU = 0.0;
		
		graph.forEdgesOf(u, [&](node , node v, edgeweight w) {
			const index cv = partition[v], ncv = nextPartition[v];
			const bool previouslyCut = cu != cv, nowCut = ncu != ncv;
			if (previouslyCut) {
				cutUpdateCU -= w;
			}
			if (nowCut) {
				cutUpdateNCU += w;
			}
			
			// aggregate cut update for v, if v was not moved. if v was moved, it will aggregate the updates itself
			if (cv == ncv && previouslyCut != nowCut) {
				if (!cutUpdates.indexIsUsed(cv)) {
					cutUpdates.insert(cv, 0.0);
				}
				const double d = nowCut ? w : -w;
				cutUpdates[cv] += d;
				otherCuts += d;
			}
			
			volU += w;
			if (u == v) {
				volU += w;
			}
		});
	
		
		volumeUpdates[cu] -= volU;
		volumeUpdates[ncu] += volU;
		
		if (!cutUpdates.indexIsUsed(cu)) {
			cutUpdates.insert(cu, 0.0);
		}
		cutUpdates[cu] += cutUpdateCU;
		
		if (!cutUpdates.indexIsUsed(ncu)) {
			cutUpdates.insert(ncu, 0.0);
		}
		cutUpdates[ncu] += cutUpdateNCU;
		
		totalCutUpdate += cutUpdateCU + cutUpdateNCU + otherCuts;
	}

	
	// apply aggregated updates to cut and volume
	for (const index& c : cutUpdates.insertedIndexes()) {
		assert(cutUpdates.indexIsUsed(c));
		
		if (volumeUpdates[c] != 0.0) {
			// every cluster is updated at most once by every thread. if it becomes a bottleneck, we can still recursively merge clearlists in parallel
			#pragma omp atomic update
				clusterVolume[c] += volumeUpdates[c];
			volumeUpdates[c] = 0.0;
		}
		
		if (cutUpdates[c] != 0.0) {
			#pragma omp atomic update
				clusterCut[c] += cutUpdates[c];
		}
		cutUpdates.resetEntry(c);
	}
	cutUpdates.clearIndexes();

	#pragma omp atomic update
		totalCut += totalCutUpdate;
}


void LouvainMapEquation::calculateInitialClusterCutAndVolume() {
	totalCut = additionalCut; totalVolume = additionalVolume;
	#pragma omp parallel if (graph.upperNodeIdBound() > 1e6)
	{
		double tCut = 0, tVol = 0;
		#pragma omp for schedule (guided)
		for (node u = 0; u < graph.upperNodeIdBound(); ++u) {
			if (graph.hasNode(u)) {
				graph.forEdgesOf(u, [&](node , node v, edgeweight ew) {
					if (u != v) {
						clusterCut[u] += ew;
					}
					else {
						ew *= 2;
					}
					clusterVolume[u] += ew;
				});
			}
			tCut += clusterCut[u];
			tVol += clusterVolume[u];
		}
		
		#pragma omp atomic update
			totalCut += tCut;
		#pragma omp atomic update
			totalVolume += tVol;
	}
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
	// free some memory
	clusterVolume.clear();
	clusterVolume.shrink_to_fit();
	clusterCut.clear();
	clusterCut.shrink_to_fit();

	ParallelPartitionCoarsening coarsening(graph, partition, true, graph.numberOfNodes() > 1e6);
	coarsening.run();
	const Graph& metaGraph = coarsening.getCoarseGraph();
	const auto& fineToCoarseMapping = coarsening.getFineToCoarseNodeMapping();

	INFO("Run hierarchical with ", metaGraph.numberOfNodes(), " clusters (from ", graph.numberOfNodes(), " nodes)");

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

void LouvainMapEquation::checkUpdatedCutsAndVolumesAgainstRecomputation() {
	std::vector<double> cut(graph.upperNodeIdBound(), 0.0), vol(graph.upperNodeIdBound(), 0.0);
	double tCut = additionalCut, tVol = additionalVolume;
	graph.forNodes([&](node u) {
		double volU = 0.0, cutU = 0.0;
		const index cu = partition[u];
		graph.forEdgesOf(u, [&](node, node v, edgeweight weight) {
			if (cu != partition[v]) {
				cutU += weight;
			}
			
			if (u == v)
				weight *= 2;
			volU += weight;
		});
		
		vol[partition[u]] += volU;
		tVol += volU;
		
		cut[partition[u]] += cutU;
		tCut += cutU;
	});
	
	graph.forNodes([&](node u) {
		const index cu = partition[u];
		assert(cu == nextPartition[u]);
		assert(vol[cu] == clusterVolume[cu]);
		assert(cut[cu] == clusterCut[cu]);
	});
	assert(tVol == totalVolume);
	assert(tCut == totalCut);
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

