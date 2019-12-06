/*
 * LPPotts.cpp
 *
 * Created on: 2019-01-14
 * Author: Armin Wiebigke
 */

#include <omp.h>
#include <algorithm>

#include "LPPotts.h"
#include "../Globals.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Timer.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/UniformRandomSelector.h"

namespace NetworKit {

LPPotts::LPPotts(const Graph &G, double alpha, count theta, count maxIterations,
				 bool parallelPropagation)
		: CommunityDetectionAlgorithm(G), alpha(alpha), updateThreshold(theta),
		  maxIterations(maxIterations), parallelPropagation(parallelPropagation) {
}

LPPotts::LPPotts(const Graph &G, const Partition &baseClustering, double alpha, count theta,
                 count maxIterations, bool parallelPropagation)
		: CommunityDetectionAlgorithm(G, baseClustering), alpha(alpha), updateThreshold(theta),
		  maxIterations(maxIterations), parallelPropagation(parallelPropagation) {

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
	}

	count nUpdated; // number of nodes which have been updated in last iteration
	nUpdated = n; // all nodes have new labels -> first loop iteration runs
	nIterations = 0; // number of iterations
	std::vector<bool> activeNodes(z, true); // record if node must be processed

	Partition secondPartition = result;
	Partition *nextResult = &result;
	if (parallelPropagation)
		nextResult = &secondPartition;

	std::vector<count> globalLabelCounts(z, 1), glc_2(z,
													  1); // record the number of nodes for each label
	std::vector<count> *nextGlobalLabelCounts = &globalLabelCounts;
	if (parallelPropagation)
		nextGlobalLabelCounts = &glc_2;

	Aux::UniformRandomSelector selector;
	Aux::Timer runtime;

	// propagate labels
	while ((nUpdated > this->updateThreshold) && (nIterations <
												  maxIterations)) { // as long as a label has changed... or maximum iterations reached
		runtime.start();
		nIterations += 1;
		DEBUG("[BEGIN] LabelPropagation: iteration #", nIterations);

		// reset updated
		nUpdated = 0;

		G.forNodesInRandomOrder([&](node u) {
			if ((activeNodes[u]) && (G.degree(u) > 0)) {
				using label = index; // a label is the same as a cluster id
				std::map<label, count> neighborLabelCounts; // neighborLabelCounts maps label -> frequency in the neighbors

				// Count the labels of the neighbors
				G.forNeighborsOf(u, [&](node w) {
					label lw = result.subsetOf(w);
					++neighborLabelCounts[lw];
				});

				// Evaluate each label
				std::map<label, double> labelWeights;
				for (const auto &labelCount : neighborLabelCounts) {
					label l = labelCount.first;
					count local = labelCount.second;
					labelWeights[l] = local - alpha * (globalLabelCounts[l] - local);
				}

				// Get best label
				selector.reset();
				label bestLabel = none;
				double bestWeight = -std::numeric_limits<double>::max();
				for (auto current : labelWeights) {
					if (current.second > bestWeight) {
						bestWeight = current.second;
						bestLabel = current.first;
						selector.reset();
					} else if (current.second == bestWeight) {
						if (selector.addElement())
							bestLabel = current.first;
					}
				}
				assert(bestLabel != none);

				label curLabel = result.subsetOf(u);
				if (curLabel != bestLabel) { // UPDATE
					--(*nextGlobalLabelCounts)[curLabel];
					++(*nextGlobalLabelCounts)[bestLabel];
					nextResult->moveToSubset(bestLabel, u);
					++nUpdated;
					G.forNeighborsOf(u, [&](node w) {
						activeNodes[w] = true;
					});
				} else {
					activeNodes[u] = false;
				}

			} else {
				// node is isolated
			}
		});

		if (parallelPropagation) {
			result = *nextResult;
			globalLabelCounts = *nextGlobalLabelCounts;
		}

		// for each while loop iteration...

		runtime.stop();
		this->timing.push_back(runtime.elapsedMilliseconds());
		DEBUG("[DONE] LabelPropagation: iteration #", nIterations,
			  " - updated ", nUpdated, " labels, time spent: ",
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


LPPottsFactory::LPPottsFactory(double alpha, count theta, count maxIterations, bool parallelPropagation)
		: alpha(alpha), theta(theta), maxIterations(maxIterations), parallelPropagation(parallelPropagation) {
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
