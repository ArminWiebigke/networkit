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
	
	// chunkBorders = fixed number of nodes, or degree sum ?
	const size_t chunkSize = std::min(static_cast<size_t>(10000 * Aux::getCurrentNumberOfThreads()), std::max(1UL, nodes.size() / 5));
	const size_t numberOfChunks = 1 + nodes.size() / chunkSize;
	std::vector<size_t> chunkBorders = Aux::Parallel::Chunking::getChunkBorders(nodes.size(), numberOfChunks);
	
	std::vector< SparseVector<double> > ets_neighborClusterWeights(Aux::getMaxNumberOfThreads(), SparseVector<double>(graph.upperNodeIdBound(), 0.0));
	std::vector< std::vector<bool> > ets_isNodeInCurrentChunk(Aux::getMaxNumberOfThreads(), std::vector<bool>(graph.upperNodeIdBound(), false));
	std::vector< NeighborCaches > ets_neighborCaches(Aux::getMaxNumberOfThreads(), NeighborCaches(chunkSize + 1, NeighborCache(10)));
	
	for (count iteration = 0; iteration < maxIterations; ++iteration) {
		handler.assureRunning();

		DEBUG("Iteration ", iteration);
#ifndef NDEBUG
		DEBUG("Map equation is ", mapEquation());
#endif

		timer.start();
		std::shuffle(nodes.begin(), nodes.end(), Aux::Random::getURNG());
		timer.stop();
		INFO("shuffle ", timer.elapsedMilliseconds(), " ms");

		count numberOfNodesMoved = 0;
		
		timer.start();
		#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			
			SparseVector<double>& neighborClusterWeights = ets_neighborClusterWeights[tid];
			//neighborClusterWeights.resize(graph.upperNodeIdBound(), 0.0);
			
			std::vector<Move> moves;
			NeighborCaches& neighborCaches = ets_neighborCaches[tid];
			index numUsedCaches = 0;
			
			std::vector<bool>& isNodeInCurrentChunk = ets_isNodeInCurrentChunk[tid];
			
			for (index i = 0; i < chunkBorders.size() - 1; ++i) {
				const index firstInvalid = chunkBorders[i + 1];
				
				for (index j = chunkBorders[i]; j < firstInvalid; ++j) {
					isNodeInCurrentChunk[ nodes[j] ] = true;
				}
				
				// find moves
				numUsedCaches = 0;
				#pragma omp for
				for (index j = chunkBorders[i]; j < firstInvalid; ++j) {
					const node u = nodes[j];
					if (numUsedCaches == neighborCaches.size()) {
						neighborCaches.emplace_back(10);	// construct neighbor cache with space for 10 entries
					}
					tryLocalMove(u, numUsedCaches /* by reference, updated in function */, neighborCaches[numUsedCaches], neighborClusterWeights, moves, isNodeInCurrentChunk);
				}
				// implicit barrier at the end of the loop
				
				aggregateAndApplyCutAndVolumeUpdates(moves, neighborCaches);
				#pragma omp atomic update
				numberOfNodesMoved += moves.size();
				moves.clear();
				
				for (index j = chunkBorders[i]; j < firstInvalid; ++j) {
					isNodeInCurrentChunk[ nodes[j] ] = false;
				}
				#pragma omp barrier
				
				#ifndef NDEBUG
					#pragma omp single
					{
						checkUpdatedCutsAndVolumesAgainstRecomputation();
					}
					// don't start the next round of moves, before the old ones were validated
					#pragma omp barrier
				#endif
				
			}
		}
		timer.stop();
		INFO("Move iteration ", iteration, " took ", timer.elapsedMilliseconds(), " ms. Moved ", numberOfNodesMoved, " nodes");
		
		clusteringChanged |= numberOfNodesMoved > 0;
		if (numberOfNodesMoved == 0) {
			break;
		}
	}

	handler.assureRunning();
	if (hierarchical && clusteringChanged) {
		runHierarchical();
	}
	hasRun = true;
}

void LouvainMapEquation::aggregateAndApplyCutAndVolumeUpdates(std::vector<Move>& moves, NeighborCaches& neighborCaches) {
	// can get rid of atomic updates by aggregating everything in clearlists
	
	for (Move& move : moves) {
		const index originCluster = move.originCluster;
		const index targetCluster = move.targetCluster;
		const node u = move.movedNode;
		
		// apply volume updates
		#pragma omp atomic update
			clusterVolume[originCluster] -= move.volume;
		#pragma omp atomic update
			clusterVolume[targetCluster] += move.volume;
			
		double 	originClusterCutUpdate = move.cutUpdateToOriginCluster,
				targetClusterCutUpdate = move.cutUpdateToTargetCluster;
			
		// correct the cut for potentially moved neighbors
		// the already applied cut updates assumed none of them were moved
		for (NeighborInChunk& cn : neighborCaches[move.cacheID]) {
			const node neighbor = cn.neighbor;
			assert(u < neighbor);
			const index originClusterOfNeighbor = cn.oldCluster;
			const index targetClusterOfNeighbor = nextPartition[neighbor];
			if (targetClusterOfNeighbor != originClusterOfNeighbor) {
				const double w2 = 2 * cn.weightToNeighbor;
			
				if (originCluster == originClusterOfNeighbor) {
					originClusterCutUpdate -= w2;
				} else if (targetClusterOfNeighbor == originCluster) {
					originClusterCutUpdate += w2;
				}
				
				if (targetClusterOfNeighbor == targetCluster) {
					targetClusterCutUpdate -= w2;
				} else if (originClusterOfNeighbor == targetCluster) {
					targetClusterCutUpdate += w2;
				}
			}
		}
		
		// apply cut updates
		#pragma omp atomic update
		clusterCut[originCluster] += originClusterCutUpdate;
		#pragma omp atomic update
		clusterCut[targetCluster] += targetClusterCutUpdate;
		#pragma omp atomic update
		totalCut += originClusterCutUpdate + targetClusterCutUpdate;
		
		
		// write new partition
		partition[u] = targetCluster;
	}
}





// for every node. store its neighbors that are in the current chunk, and their old cluster IDs, and edge weights
bool LouvainMapEquation::tryLocalMove(node u, index& cacheID, std::vector<NeighborInChunk>& cachedNeighbors, SparseVector<double>& neighborClusterWeights,
									  std::vector<Move>& moves, std::vector<bool>& isNodeInCurrentChunk) {
	// Find neighbor clusters
	double vol = 0;
	double loop = 0;
	double weightToCurrent = 0;
	const index currentCluster = partition[u];
	cachedNeighbors.clear();
	
	graph.forEdgesOf(u, [&](node, node v, edgeweight weight) {
		vol += weight;
		if (u != v) {
			const index neighborCluster = partition[v];
			
			if (u < v && isNodeInCurrentChunk[v]) {	// TODO global bitvector or local? global has to be synchronized with an extra barrier; got anything better?
				cachedNeighbors.emplace_back(v, neighborCluster, weight);
			}
			
			if (neighborCluster == currentCluster) {
				weightToCurrent += weight;
			} else {
				if (!neighborClusterWeights.indexIsUsed(neighborCluster))
					neighborClusterWeights.insert(neighborCluster, 0);
				neighborClusterWeights[neighborCluster] += weight;
			}
		} else {
			loop += weight;
			vol += weight;
		}
	});
	
	assert(vol == graph.weightedDegree(u));

	if (neighborClusterWeights.size() > 0) {
		// Calculate best cluster
		index targetCluster = currentCluster;
		double weightToTargetCluster = weightToCurrent;
		double bestChange = fitnessChange(u, vol, loop, currentCluster, currentCluster, weightToCurrent, weightToCurrent);
		for (index neighborCluster : neighborClusterWeights.insertedIndexes()) {
			const double neighborClusterWeight = neighborClusterWeights[neighborCluster];
			const double change = fitnessChange(u, vol, loop, currentCluster, neighborCluster, neighborClusterWeight, weightToCurrent);
			if (change < bestChange || (change == bestChange && neighborCluster < targetCluster && targetCluster != currentCluster)) {
				bestChange = change;
				targetCluster = neighborCluster;
				weightToTargetCluster = neighborClusterWeight;
			}
		}
		
		neighborClusterWeights.reset();
		
		// don't apply move yet. just save it.
		if (targetCluster != currentCluster) {
		 	const double cutUpdateToCurrentCluster = -vol + 2 * weightToCurrent + 2 * loop;
		 	const double cutUpdateToTargetCluster = vol - 2 * weightToTargetCluster - 2 * loop;
			moves.emplace_back(u, vol, cacheID++, currentCluster, targetCluster, cutUpdateToCurrentCluster, cutUpdateToTargetCluster);
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

	assert(metaGraph.numberOfNodes() < graph.numberOfNodes());
	
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
LouvainMapEquation::fitnessChange(node, double degree, double loopWeight, index currentCluster, index targetCluster,
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
		if (x > 0.0) {
			x /= totalVolume;
			x *= std::log(x);
		} else {
			x = 0.0;
		}
	};

	// for unit edge weights: wouldn't it make sense to precompute a lookup table for these logarithms?
	
	normalizeAndPLogP(totalCutNew);
	normalizeAndPLogP(targetClusterCutNew);
	normalizeAndPLogP(targetClusterCutCurrent);
	normalizeAndPLogP(targetCutPlusVolumeNew);
	normalizeAndPLogP(targetCutPlusVolumeCurrent);

	return totalCutNew + ((targetCutPlusVolumeNew - targetCutPlusVolumeCurrent) - (2 * (targetClusterCutNew - targetClusterCutCurrent)));
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

