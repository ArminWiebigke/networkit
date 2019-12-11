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

#include <networkit/community/LouvainMapEquation.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Timer.hpp>

namespace NetworKit {

LouvainMapEquation::LouvainMapEquation(const Graph &graph, bool hierarchical, count maxIterations,
									   bool parallel, ParallelizationType parallelizationType,
									   double additionalCut, double additionalVolume)
		: parallel(parallel), parallelizationType(parallelizationType),
		  graph(graph), hierarchical(hierarchical), maxIterations(maxIterations),
		  partition(graph.upperNodeIdBound()),
		  clusterCut(graph.upperNodeIdBound()),
		  clusterVolume(graph.upperNodeIdBound()),
		  additionalCut(additionalCut),
		  additionalVolume(additionalVolume),
		  totalCut(additionalCut),
		  totalVolume(additionalVolume),
		  locks(parallel && parallelizationType == ParallelizationType::RelaxMap ? graph.upperNodeIdBound() : 0),
		  nextPartition(parallel && parallelizationType == ParallelizationType::SynchronousLocalMoving ? graph.upperNodeIdBound() : 0),
		  // parallel variants initialize the inner vectors in parallel
		  ets_neighborClusterWeights(parallel ? Aux::getMaxNumberOfThreads() : 1, SparseVector<double>(parallel ? 0 : graph.upperNodeIdBound(), 0.0)),
		  ets_isNodeInCurrentChunk(parallel && parallelizationType == ParallelizationType::SynchronousLocalMoving ? Aux::getMaxNumberOfThreads() : 0),
		  ets_neighborCaches(parallel && parallelizationType == ParallelizationType::SynchronousLocalMoving? Aux::getMaxNumberOfThreads() : 0)
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
				if (parallel && parallelizationType == ParallelizationType::SynchronousLocalMoving) {
					nextPartition.remove(u);
				}
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
	
	count numberOfNodesMoved = 1;
	for (count iteration = 0; iteration < maxIterations && numberOfNodesMoved > 0; ++iteration) {
		handler.assureRunning();

		DEBUG("Iteration ", iteration);
#ifndef NDEBUG
		DEBUG("Map equation is ", mapEquation());
#endif

		timer.start();
		std::shuffle(nodes.begin(), nodes.end(), Aux::Random::getURNG());
		timer.stop();
		INFO("shuffle ", timer.elapsedMilliseconds(), " ms");
		
		timer.start();
		if (parallel && parallelizationType == ParallelizationType::SynchronousLocalMoving) {
			numberOfNodesMoved = synchronousLocalMoving(nodes, iteration);
		} else {
			numberOfNodesMoved = localMoving(nodes, iteration);
		}
		timer.stop();
		INFO("Move iteration ", iteration, " took ", timer.elapsedMilliseconds(), " ms. Moved ", numberOfNodesMoved, " nodes");
		
		clusteringChanged |= numberOfNodesMoved > 0;
	}

	handler.assureRunning();
	if (hierarchical && clusteringChanged) {
		runHierarchical();
	}
	hasRun = true;
}

count LouvainMapEquation::localMoving(std::vector<node>& nodes, count iteration) {
	
	// dummies, since SLM implementation needs more datastructures
	index dummyCacheID;
	std::vector<NeighborInChunk> dummyCachedNeighbors;
	std::vector<Move> dummyMoves;
	std::vector<bool> dummyIsNodeInCurrentChunk;
	
	count nodesMoved = 0;
	if (parallel) {
		#pragma omp parallel
		{
			count nm = 0;
			int tid = omp_get_thread_num();
			SparseVector<double>& neighborClusterWeights = ets_neighborClusterWeights[tid];

			if (iteration == 0) {
				neighborClusterWeights = SparseVector<double>(graph.upperNodeIdBound(), 0.0);
			}
			
			#pragma omp for schedule(guided) nowait
			for (index i = 0; i < nodes.size(); ++i) {
				if (tryLocalMove<true, false>(nodes[i], neighborClusterWeights, /* dummies */ dummyCacheID, dummyCachedNeighbors, dummyMoves, dummyIsNodeInCurrentChunk)) {
					nm += 1;
				}
			}

			#pragma omp atomic update
			nodesMoved += nm;
		}
	} else {
		SparseVector<double>& neighborClusterWeights = ets_neighborClusterWeights[0];
		for (node u : nodes) {
			if (tryLocalMove<false, false>(u, neighborClusterWeights, /* dummies */ dummyCacheID, dummyCachedNeighbors, dummyMoves, dummyIsNodeInCurrentChunk)) {
				nodesMoved += 1;
			}
		}
	}
	return nodesMoved;
}


count LouvainMapEquation::synchronousLocalMoving(std::vector<NetworKit::node>& nodes, count iteration) {
	// Estimate the number of nodes that will be moved in this iteration and make chunkSize dependent on that! Or implement active nodesets
	const size_t chunkSize = std::min(static_cast<size_t>(10000 * Aux::getCurrentNumberOfThreads()), std::max(1UL, nodes.size() / 5));
	const size_t numberOfChunks = 1 + nodes.size() / chunkSize;
	std::vector<size_t> chunkBorders = Aux::Parallel::Chunking::getChunkBorders(nodes.size(), numberOfChunks);
	
	count numberOfNodesMoved = 0;
	
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		
		SparseVector<double>& neighborClusterWeights = ets_neighborClusterWeights[tid];
		
		std::vector<Move> moves;
		NeighborCaches& neighborCaches = ets_neighborCaches[tid];
		index numUsedCaches = 0;
		
		std::vector<bool>& isNodeInCurrentChunk = ets_isNodeInCurrentChunk[tid];
		
		if (iteration == 0) {	// with thread pinning, this ensures the datastructures are allocated on the right socket
			isNodeInCurrentChunk = std::vector<bool>(graph.upperNodeIdBound(), false);
			neighborClusterWeights = SparseVector<double>(graph.upperNodeIdBound(), 0.0);
			neighborCaches = NeighborCaches(chunkSize + 1 , NeighborCache(20));
		}
		
		
		for (index i = 0; i < chunkBorders.size() - 1; ++i) {
			const index firstInvalid = chunkBorders[i + 1];
			
			for (index j = chunkBorders[i]; j < firstInvalid; ++j) {
				isNodeInCurrentChunk[ nodes[j] ] = true;
			}
			
			// find moves
			numUsedCaches = 0;
			#pragma omp for schedule (guided)
			for (index j = chunkBorders[i]; j < firstInvalid; ++j) {
				const node u = nodes[j];
				if (numUsedCaches == neighborCaches.size()) {
					neighborCaches.emplace_back(10);	// construct neighbor cache with space for 10 entries
				}
				tryLocalMove<true, true>(u, neighborClusterWeights, numUsedCaches /* by reference, updated in function */, neighborCaches[numUsedCaches], moves, isNodeInCurrentChunk);
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
	return numberOfNodesMoved;
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

template<bool pparallel, bool synchronous>
bool LouvainMapEquation::tryLocalMove(node u, SparseVector<double>& neighborClusterWeights,
												 /* SLM specifics */
												 index& cacheID, NeighborCache& cachedNeighbors, std::vector<Move>& moves, std::vector<bool>& isNodeInCurrentChunk) {
	// Find neighbor clusters
	double vol = 0;
	double loop = 0;
	double weightToCurrent = 0;
	const index currentCluster = partition[u];
	
	if /* constexpr */ (synchronous) {
		cachedNeighbors.clear();
	}
	
	graph.forEdgesOf(u, [&](node, node v, edgeweight weight) {
		vol += weight;
		if (u != v) {
			const index neighborCluster = partition[v];
			
			if /* constexpr */ (synchronous) {
				if (u < v && isNodeInCurrentChunk[v]) {
					cachedNeighbors.emplace_back(v, neighborCluster, weight);
				}
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
	
	assert(vol == graph.weightedDegree(u, true));

	if (neighborClusterWeights.size() > 0) {
		// Calculate best cluster
		index targetCluster = currentCluster;
		double weightToTargetCluster = weightToCurrent;
		double bestChange = fitnessChange(u, vol, loop, currentCluster, currentCluster, weightToCurrent, weightToCurrent, totalCut);
		for (index neighborCluster : neighborClusterWeights.insertedIndexes()) {
			const double neighborClusterWeight = neighborClusterWeights[neighborCluster];
			const double change = fitnessChange(u, vol, loop, currentCluster, neighborCluster, neighborClusterWeight, weightToCurrent, totalCut);
			if (change < bestChange || (change == bestChange && neighborCluster < targetCluster && targetCluster != currentCluster)) {
				bestChange = change;
				targetCluster = neighborCluster;
				weightToTargetCluster = neighborClusterWeight;
			}
		}
		
		neighborClusterWeights.reset();
		
		if (targetCluster != currentCluster) {
			if /* constexpr */ (synchronous) {
				assert(parallel);
				// save move
				const double cutUpdateToCurrentCluster = -vol + 2 * weightToCurrent + 2 * loop;
				const double cutUpdateToTargetCluster = vol - 2 * weightToTargetCluster - 2 * loop;
				moves.emplace_back(u, vol, cacheID++, currentCluster, targetCluster, cutUpdateToCurrentCluster, cutUpdateToTargetCluster);
				nextPartition[u] = targetCluster;
				return true;
			} else {
				// perform move directly
				return performMove<pparallel>(u, vol, loop, currentCluster, targetCluster, weightToTargetCluster, weightToCurrent);
			}
		}
	}
	return false;
}

template<bool parallel>
bool LouvainMapEquation::performMove(node u, double degree, double loopWeight, node currentCluster, node targetCluster, double weightToTarget, double weightToCurrent) {
	bool moved = true;
	if /* constexpr */ (parallel) {
		assert(parallelizationType == ParallelizationType::RelaxMap);
		
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
		
		const double totalCutCurrently = totalCut;
		const double fitnessCurrent = fitnessChange(u, degree, loopWeight, currentCluster, currentCluster, weightToCurrent, weightToCurrent, totalCutCurrently);
		const double fitnessTarget = fitnessChange(u, degree, loopWeight, currentCluster, targetCluster, weightToTarget, weightToCurrent, totalCutCurrently);
		if (fitnessTarget >= fitnessCurrent) {
			moved = false;
		}
	}
	
	if (moved) {
		double cutDifferenceCurrent = 2 * weightToCurrent - degree + 2 * loopWeight;
		double cutDifferenceTarget = degree - 2 * weightToTarget - 2 * loopWeight;
		clusterCut[currentCluster] += cutDifferenceCurrent;
		clusterCut[targetCluster] += cutDifferenceTarget;
		clusterVolume[currentCluster] -= degree;
		clusterVolume[targetCluster] += degree;
		partition.moveToSubset(targetCluster, u);
		
		if /*constexpr */ (parallel) {
			totalCut += cutDifferenceCurrent + cutDifferenceTarget;
		} else {
			#pragma omp atomic
			totalCut += cutDifferenceCurrent + cutDifferenceTarget;
		}
	}
	
	if /* constexpr */ (parallel) {
		// unlock clusters again
		locks[ std::max(currentCluster, targetCluster) ].unlock();
		locks[ std::min(currentCluster, targetCluster) ].unlock();
	}
	
	return moved;
}

void LouvainMapEquation::runHierarchical() {
	assert(partition.numberOfSubsets() < partition.numberOfElements());
	// free some memory
	clusterVolume.clear();
	clusterVolume.shrink_to_fit();
	clusterCut.clear();
	clusterCut.shrink_to_fit();

	ParallelPartitionCoarsening coarsening(graph, partition, true, parallel && graph.numberOfNodes() > 1e6);
	coarsening.run();
	const Graph& metaGraph = coarsening.getCoarseGraph();
	const auto& fineToCoarseMapping = coarsening.getFineToCoarseNodeMapping();

	assert(metaGraph.numberOfNodes() < graph.numberOfNodes());
	
	INFO("Run hierarchical with ", metaGraph.numberOfNodes(), " clusters (from ", graph.numberOfNodes(), " nodes)");

	constexpr bool forceParallelism = false;	// for tests only
	LouvainMapEquation recursion(metaGraph, true, maxIterations, parallel && (forceParallelism || metaGraph.numberOfNodes() > 1e4), parallelizationType, additionalCut, additionalVolume);
	recursion.run();
	const Partition& metaPartition = recursion.getPartition();

	graph.forNodes([&](node u) {
		partition[u] = metaPartition[fineToCoarseMapping[u]];
	});
}

const Partition& LouvainMapEquation::getPartition() const {
	return partition;
}

void LouvainMapEquation::calculateInitialClusterCutAndVolume() {
	totalCut = additionalCut; totalVolume = additionalVolume;
	
	if (parallel) {
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
	} else {
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
			totalCut += clusterCut[u];
			totalVolume += clusterVolume[u];
		}
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

LouvainMapEquationFactory::LouvainMapEquationFactory(bool hierarchical, count maxIterations, std::string parallelization)
		: hierarchical(hierarchical), maxIterations(maxIterations), parallelization(std::move(parallelization)) {
}

ClusteringFunction LouvainMapEquationFactory::getFunction() const {
	bool hierarchicalCopy = hierarchical;
	count maxIterationsCopy = maxIterations;
	
	// Unknown parallelization options default to no parallelism!
	bool parallel = parallelization == "SynchronousLocalMoving" || parallelization == "RelaxMap";
	LouvainMapEquation::ParallelizationType parallelizationType = LouvainMapEquation::ParallelizationType::RelaxMap;
	if (parallelization == "SynchronousLocalMoving") {
		parallelizationType = LouvainMapEquation::ParallelizationType::SynchronousLocalMoving;
	}
	
	return [hierarchicalCopy, maxIterationsCopy, parallel, parallelizationType](const Graph &graph) {
		LouvainMapEquation algo(graph, hierarchicalCopy, maxIterationsCopy, parallel, parallelizationType);
		algo.run();
		return algo.getPartition();
	};
}

}

