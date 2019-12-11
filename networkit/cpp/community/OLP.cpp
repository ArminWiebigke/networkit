/*
 * OLP.cpp
 *
 * Created on: 2019-01-28
 * Author: Armin Wiebigke
 */

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/auxiliary/Random.hpp>
#include "networkit/auxiliary/UniformRandomSelector.hpp"
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/community/OLP.hpp>

namespace NetworKit {

OLP::OLP(const Graph &G, count k, count theta, count maxIterations)
		: G(G), maxLabels(k), updateThreshold(theta),
		  maxIterations(maxIterations) {
}

void OLP::run() {
	// set unique label for each node if no baseClustering was given
	index z = G.upperNodeIdBound();
	if (hasRun || result.numberOfElements() != z) {
		result = Partition(z);
		result.allToSingletons();
	}

	index labelBound = result.upperBound();

	count n = G.numberOfNodes();
	// update threshold heuristic
	if (updateThreshold == none) {
		updateThreshold = (count) (n / 1e5);
	}

	count nUpdated; // number of nodes which have been updated in last iteration
	nUpdated = n; // all nodes have new labels -> first loop iteration runs
	nIterations = 0; // number of iterations
	std::vector<bool> activeNodes(z, true); // record if node must be processed

	std::vector<count> globalLabelCounts(labelBound,
										 1); // record the number of nodes for each label

	Aux::UniformRandomSelector selector;
	Aux::Timer runtime;
	timing.clear();
	Aux::SignalHandler handler;

	// propagate labels
	while ((nUpdated > this->updateThreshold) && (nIterations <
												  maxIterations)) { // as long as a label has changed... or maximum iterations reached
		handler.assureRunning();
		runtime.start();
		nIterations += 1;
		DEBUG("[BEGIN] LabelPropagation: iteration #", nIterations);

		// reset updated
		nUpdated = 0;

		G.forNodesInRandomOrder([&](node u) {
			if ((activeNodes[u]) && (G.degree(u) > 0)) {
				using label = index; // a label is the same as a cluster id
				std::map<label, double> neighborLabelCounts; // neighborLabelCounts maps label -> (weighted) frequency in the neighbors

				// Count the labels of the neighbors
				G.forNeighborsOf(u, [&](node w, edgeweight weight) {
					auto labels = result.subsetsOf(w);
					for (label l : labels)
						neighborLabelCounts[l] += weight;
				});

				std::set<label> bestLabels;

				// Get best labels
				if (neighborLabelCounts.size() <= maxLabels) {
					for (auto &label_pair : neighborLabelCounts) {
						bestLabels.insert(label_pair.first);
					}
				} else {
					// Sort the labels by count
					std::vector<std::pair<label, count>> pairs;
					for (auto &label_pair : neighborLabelCounts)
						pairs.emplace_back(label_pair);
					sort(pairs.begin(), pairs.end(), [](std::pair<label, count> &a,
														std::pair<label, count> &b) {
							 return a.second > b.second;
						 }
					);
					// Shuffle labels with same count
					if (pairs[maxLabels - 1].second == pairs[maxLabels].second) {
						count sameCount = pairs[maxLabels].second;
						int i = 0;
						while (pairs[i].second > sameCount) ++i;
						int j = i;
						while (j < pairs.size() && pairs[j].second == sameCount) ++j;
						std::shuffle(pairs.begin() + i, pairs.begin() + j,
									 Aux::Random::getURNG());
					}
					for (std::size_t i = 0; i < maxLabels; ++i) {
						bestLabels.insert(pairs[i].first);
					}
				}

				assert(bestLabels.size() <= maxLabels);

				auto curLabels = result.subsetsOf(u);
				assert(curLabels.size() <= maxLabels);
				if (curLabels != bestLabels) { // UPDATE
					for (label l : curLabels) {
						assert(l < labelBound);
						--globalLabelCounts[l];
						result.removeFromSubset(l, u);
					}
					for (label l : bestLabels) {
						assert(l < labelBound);
						++globalLabelCounts[l];
						result.addToSubset(l, u);
					}
					assert(result.subsetsOf(u).size() <= maxLabels);
					++nUpdated;
					G.forNeighborsOf(u, [&](node v) {
						activeNodes[v] = true;
					});
				} else {
					activeNodes[u] = false;
				}

			} else {
				// node is isolated
			}
		});

		// for each while loop iteration...

		runtime.stop();
		this->timing.push_back(runtime.elapsedMilliseconds());
		DEBUG("[DONE] LabelPropagation: iteration #", nIterations,
			  " - updated ", nUpdated, " labels, time spent: ",
			  runtime.elapsedTag());
	} // end while

	// Discard communities of size 4 or less
	count min_size = 5;
	std::vector<std::vector<node>> communities(result.upperBound());
	G.forNodes([&](node u) {
		for (index c : result.subsetsOf(u)) {
			if (communities[c].size() < min_size)
				communities[c].push_back(u);
		}
	});
	for (index c = 0; c < communities.size(); ++c) {
		if (communities[c].size() < min_size) {
			for (node u : communities[c])
				result.removeFromSubset(c, u);
		}
	}

	hasRun = true;
}

std::string OLP::toString() const {
	std::stringstream strm;
	strm << "OLP";
	return strm.str();
}


void OLP::setUpdateThreshold(count th) {
	this->updateThreshold = th;
}


count OLP::numberOfIterations() {
	return this->nIterations;
}


std::vector<count> OLP::getTiming() {
	return this->timing;
}

Cover OLP::getCover() {
	return result;
}

} /* namespace NetworKit */
