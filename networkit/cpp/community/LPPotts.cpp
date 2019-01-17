/*
 * LPPotts.cpp
 *
 *  Created on: 07.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "LPPotts.h"

#include <omp.h>
#include "../Globals.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Timer.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

LPPotts::LPPotts(const Graph &G, double alpha, count theta, count maxIterations)
        : CommunityDetectionAlgorithm(G), alpha(alpha), updateThreshold(theta),
          maxIterations(maxIterations) {
}

//LPPotts::LPPotts(const Graph &G, const Partition baseClustering, count theta)
//        : CommunityDetectionAlgorithm(G, baseClustering),
//          updateThreshold(theta) {
//}

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
    std::vector<count> globalLabelCounts(z, 1); // record the number of nodes for each label

    Aux::Timer runtime;

    // propagate labels
    while ((nUpdated > this->updateThreshold) && (nIterations < maxIterations)) { // as long as a label has changed... or maximum iterations reached
        runtime.start();
        nIterations += 1;
        DEBUG("[BEGIN] LabelPropagation: iteration #", nIterations);

        // reset updated
        nUpdated = 0;

        // TODO: random permumation of nodes

        G.forNodes([&](node v) {
            if ((activeNodes[v]) && (G.degree(v) > 0)) {
                using label = index; // a label is the same as a cluster id
                std::map<label, count> neighborLabelCounts; // neighborLabelCounts maps label -> frequency in the neighbors

                // Count the labels of the neighbors
                G.forNeighborsOf(v, [&](node w) {
                    label lw = result.subsetOf(w);
                    ++neighborLabelCounts[lw];
                });

                // Evaluate each label
                std::map<label, double> labelWeights;
                for (const auto &labelCount : neighborLabelCounts) {
                    label l = labelCount.first;
                    count local = labelCount.second;
                    labelWeights[l] =
                            local - alpha * (globalLabelCounts[l] - local);
                }

                // Get best label
                auto labelComparator = [](
                        const std::pair<label, double> &p1,
                        const std::pair<label, double> &p2) {
                    return p1.second < p2.second;
                };
                label best = std::max_element(labelWeights.begin(),
                                              labelWeights.end(),
                                              labelComparator
                )->first;

                if (result.subsetOf(v) != best) { // UPDATE
                    --globalLabelCounts[result.subsetOf(v)];
                    result.moveToSubset(best, v);
                    ++globalLabelCounts[best];
                    ++nUpdated;
                    G.forNeighborsOf(v, [&](node u) {
                        activeNodes[u] = true;
                    });
                } else {
                    activeNodes[v] = false;
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

} /* namespace NetworKit */
