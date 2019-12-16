/*
 * LPPotts.cpp
 *
 * Created on: 2019-01-14
 * Author: Armin Wiebigke
 */

#include <omp.h>
#include <algorithm>
#include <cmath>

#include <networkit/community/LPPotts.hpp>
#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/UniformRandomSelector.hpp>
#include <networkit/auxiliary/SparseVector.hpp>
#include <networkit/auxiliary/Parallelism.hpp>

namespace NetworKit {

LPPotts::LPPotts(const Graph &G, double alpha, count theta, count maxIterations,
                 bool parallelPropagation)
		: LPPotts(G, Partition(), alpha, theta, maxIterations, parallelPropagation) {
}

LPPotts::LPPotts(const Graph &G, const Partition &baseClustering, double alpha, count theta,
                 count maxIterations, bool parallelPropagation)
		: CommunityDetectionAlgorithm(G, baseClustering), alpha(alpha), updateThreshold(theta),
		  maxIterations(maxIterations), parallel(parallelPropagation) {
	count numThreads = 1;

	if (parallelPropagation) {
		numThreads = Aux::getMaxNumberOfThreads();
	}
	index upperNodeIdBound = G.upperNodeIdBound();
	neighborLabelCountsPerThread.resize(numThreads,
	                                    SparseVector<count>(upperNodeIdBound, none));
	labelWeightsPerThread.resize(numThreads,
	                             SparseVector<double>(upperNodeIdBound,
	                                                  -std::numeric_limits<double>::max()));
}

void LPPotts::run() {
	if (hasRun) {
		throw std::runtime_error("The algorithm has already run on the graph.");
	}

	// set unique label for each node if no baseClustering was given
	index z = G.upperNodeIdBound();
	if (result.numberOfElements() != z) {
		result = Partition(z);
		result.allToSingletons();
	}

	count n = G.numberOfNodes();
	// update threshold heuristic
	if (updateThreshold == none) {
		updateThreshold = (count) (n / 1e5);
		updateThreshold = std::max(updateThreshold, (count) 1);
	}

	count nUpdated; // number of nodes which have been updated in last iteration
	nUpdated = n; // all nodes have new labels -> first loop iteration runs
	nIterations = 0; // number of iterations

	Partition *nextResult = &result;
	Partition secondPartition;
	if (parallel) {
		secondPartition = result;
		nextResult = &secondPartition;
	}

	std::vector<count> globalLabelCounts(z,
	                                     1), globalLabelCounts2; // record the number of nodes for each label
	std::vector<count> *nextGlobalLabelCounts = &globalLabelCounts;
	if (parallel) {
		globalLabelCounts2 = globalLabelCounts;
		nextGlobalLabelCounts = &globalLabelCounts2;
	}
	std::vector<count> activeNodes(z, true); // record if node must be processed
	std::vector<count> nextActiveNodes(z, false);

	Aux::UniformRandomSelector selector;
	Aux::Timer runtime;


	// propagate labels
	while ((nUpdated > this->updateThreshold) && (nIterations < maxIterations)) {
		// as long as a label has changed... or maximum iterations reached
		runtime.start();
		nIterations += 1;
		DEBUG("[BEGIN] LabelPropagation: iteration #", nIterations);

		// reset updated
		nUpdated = 0;

#pragma omp parallel
		{

		}

		auto updateNode = [&](node u) {
			if (!activeNodes[u] || (G.degree(u) == 0)) {
				return;
			}
			index tid = parallel ? omp_get_thread_num() : 0;
			using label = index; // a label is the same as a cluster id
			// neighborLabelCounts maps label -> frequency in the neighbors
			SparseVector<count> &neighborLabelCounts = neighborLabelCountsPerThread[tid];
			assert(result.upperBound() <= z);

			// Count the labels of the neighbors
			G.forNeighborsOf(u, [&](node w) {
				label wLabel = result.subsetOf(w);
				if (!neighborLabelCounts.indexIsUsed(wLabel)) {
					neighborLabelCounts.insert(wLabel, 0);
				}
				++neighborLabelCounts[wLabel];
			});

			// Evaluate each label
			SparseVector<double> &labelWeights = labelWeightsPerThread[tid];
			for (label candidateLabel : neighborLabelCounts.insertedIndexes()) {
				count localCount = neighborLabelCounts[candidateLabel];
				double weight = localCount;
				if (alpha > 0.0) {
					weight -= alpha * (globalLabelCounts[candidateLabel] - localCount);
				}
				labelWeights.insert(candidateLabel, weight);
			}
			neighborLabelCounts.reset();

			// Get best label
			selector.reset();
			label bestLabel = none;
			double bestWeight = -std::numeric_limits<double>::max();
			for (auto neighborLabel : labelWeights.insertedIndexes()) {
				double weight = labelWeights[neighborLabel];
				if (weight > bestWeight) {
					bestWeight = weight;
					bestLabel = neighborLabel;
					selector.reset();
				} else if (weight == bestWeight) {
					if (selector.addElement())
						bestLabel = neighborLabel;
				}
			}
			assert(bestLabel != none);
			labelWeights.reset();

			// Update labels
			label curLabel = result.subsetOf(u);
			if (curLabel != bestLabel) {
				nextResult->moveToSubset(bestLabel, u);
#pragma omp critical
				{
					--(*nextGlobalLabelCounts)[curLabel];
					++(*nextGlobalLabelCounts)[bestLabel];
					++nUpdated;
					G.forNeighborsOf(u, [&](node w) {
						nextActiveNodes[w] = true;
					});
					nextActiveNodes[u] = true;
				}
			}
		};

		if (parallel) {
			G.parallelForNodes(updateNode);
			result = *nextResult;
			globalLabelCounts = *nextGlobalLabelCounts;
		} else {
			G.forNodesInRandomOrder(updateNode);
		}

		activeNodes = nextActiveNodes;
		for (count &i : nextActiveNodes)
			i = false;

		runtime.stop();
		this->timing.push_back(runtime.elapsedMilliseconds());
		DEBUG("[DONE] LabelPropagation: iteration #", nIterations, " - updated ", nUpdated,
		      " labels, time spent: ",
		      runtime.elapsedTag());
	} // end while
	hasRun = true;
}

std::string LPPotts::toString() const {
	std::stringstream strm;
	strm << "LPPotts";
	return strm.str();
}

void LPPotts::setUpdateThreshold(count th) {
	this->updateThreshold = th;
}

count LPPotts::numberOfIterations() {
	return this->nIterations;
}

std::vector<count> LPPotts::getTiming() {
	return this->timing;
}

LPPottsFactory::LPPottsFactory(double alpha, count theta, count maxIterations,
                               bool parallelPropagation)
		: alpha(alpha), theta(theta), maxIterations(maxIterations),
		  parallelPropagation(parallelPropagation) {
}

ClusteringFunction LPPottsFactory::getFunction() const {
	const double alphaCopy = alpha;
	const count thetaCopy = theta;
	const count maxIterationsCopy = maxIterations;
	const bool parallelPropagationCopy = parallelPropagation;
	return [alphaCopy, thetaCopy, maxIterationsCopy, parallelPropagationCopy](const Graph &G) {
		LPPotts lpPotts(G, alphaCopy, thetaCopy, maxIterationsCopy, parallelPropagationCopy);
		lpPotts.run();
		return lpPotts.getPartition();
	};
}
} /* namespace NetworKit */
