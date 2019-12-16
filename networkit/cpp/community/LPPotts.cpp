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
#include <networkit/community/LPDegreeOrdered.hpp>

namespace NetworKit {

LPPotts::LPPotts(const Graph &G, double alpha, count theta, count maxIterations,
                 bool parallelPropagation)
        : LPPotts(G, Partition(), alpha, theta, maxIterations, parallelPropagation) {
}

LPPotts::LPPotts(const Graph &G, const Partition &baseClustering, double alpha, count theta,
                 count maxIterations, bool parallel)
        : CommunityDetectionAlgorithm(G, baseClustering), alpha(alpha), updateThreshold(theta),
          maxIterations(maxIterations), parallel(parallel) {
    count numThreads = parallel ? Aux::getMaxNumberOfThreads() : 1;
    neighborLabelCountsPerThread.resize(numThreads,
                                        SparseVector<count>(G.upperNodeIdBound(), none));
    labelWeightsPerThread.resize(numThreads,
                                 SparseVector<double>(G.upperNodeIdBound(),
                                                      -std::numeric_limits<double>::max()));
    if (result.numberOfElements() == 0) {
        result = Partition(G.upperNodeIdBound());
        result.allToSingletons();
    }
    assert(result.numberOfElements() == G.upperNodeIdBound());

    if (updateThreshold == none) {
        updateThreshold = (count) (G.numberOfNodes() / 1e5);
        updateThreshold = std::max(updateThreshold, (count) 1);
    }
}

void LPPotts::run() {
    if (hasRun) {
        throw std::runtime_error("The algorithm has already run on the graph.");
    }

    index nodeIdBound = G.upperNodeIdBound();
    count updatedNodesCount = G.numberOfNodes();
    iteration = 0;
    Partition nextPartition;
    if (parallel) {
        nextPartition = result;
    }
    // record the number of nodes for each label
    assert(globalLabelCounts.empty());
    globalLabelCounts.resize(nodeIdBound, 1);
    std::vector<std::atomic<count>> nextGlobalLabelCountsAtomic(parallel ? nodeIdBound : 0);
    for (auto &it : nextGlobalLabelCountsAtomic) {
        std::atomic_init(&it, (count) 1);
    }
    // record if node must be processed
    std::vector<u_int8_t> activeNodes(nodeIdBound);
    G.forNodes([&](node u) { activeNodes[u] = true; });
    std::vector<u_int8_t> nextActiveNodes(nodeIdBound);
    std::vector<std::atomic<u_int8_t>> nextActiveNodesAtomic(parallel ? nodeIdBound : 0);
    for (auto &it : nextActiveNodesAtomic) {
        std::atomic_init(&it, (u_int8_t) false);
    }

    // propagate labels
    while ((updatedNodesCount > this->updateThreshold) && (iteration < maxIterations)) {
        // as long as enough labels change or maximum iterations reached
        iteration += 1;
        DEBUG("[BEGIN] LabelPropagation: iteration #", iteration);

        // reset updated
        updatedNodesCount = 0;

#pragma omp parallel if (parallel)
        {
            count localUpdatedNodes = 0;

            auto evaluateNode = [&](node u) {
                if (!activeNodes[u] || G.degree(u) == 0)
                    return;

                label bestLabel = calculateBestLabel(u);

                // Update labels
                label currentLabel = result.subsetOf(u);
                if (currentLabel != bestLabel) {
                    ++localUpdatedNodes;
                    if (parallel) {
                        nextPartition.moveToSubset(bestLabel, u);
                        nextGlobalLabelCountsAtomic[currentLabel].fetch_sub(1);
                        nextGlobalLabelCountsAtomic[bestLabel].fetch_add(1);
                        nextActiveNodesAtomic[u].store(true);
                        G.forNeighborsOf(u, [&](node w) {
                            nextActiveNodesAtomic[w].store(true);
                        });
                    } else {
                        result.moveToSubset(bestLabel, u);
                        --globalLabelCounts[currentLabel];
                        ++globalLabelCounts[bestLabel];
                        nextActiveNodes[u] = true;
                        G.forNeighborsOf(u, [&](node w) {
                            nextActiveNodes[w] = true;
                        });
                    }
                }
            };

            if (parallel) {
#pragma omp for nowait schedule(guided)
                for (node u = 0; u < nodeIdBound; ++u) {
                    evaluateNode(u);
                }
            } else {
                G.forNodesInRandomOrder(evaluateNode);
            }

#pragma omp atomic update
            updatedNodesCount += localUpdatedNodes;
        }

        if (parallel) {
            result = nextPartition;
            for (index i = 0; i < nodeIdBound; ++i) {
                activeNodes[i] = nextActiveNodesAtomic[i].load();
                nextActiveNodesAtomic[i].store(false);
                globalLabelCounts[i] = nextGlobalLabelCountsAtomic[i].load();
            }
        } else {
            for (index i = 0; i < nodeIdBound; ++i) {
                activeNodes[i] = nextActiveNodes[i];
                nextActiveNodes[i] = false;
            }
        }

//		this->timing.push_back(runtime.elapsedMilliseconds());
//		DEBUG("[DONE] LabelPropagation: iteration #", iteration, " - updated ", updatedNodesCount,
//		      " labels, time spent: ",
//		      runtime.elapsedTag());
    } // end while

    hasRun = true;
}

label LPPotts::calculateBestLabel(node u) {
    index tid = parallel ? omp_get_thread_num() : 0;

    // Count the labels of the neighbors
    SparseVector<count> &neighborLabelCounts = neighborLabelCountsPerThread[tid];
    assert(neighborLabelCounts.isClean());
    G.forNeighborsOf(u, [&](node w) {
        label neighborLabel = result.subsetOf(w);
        if (!neighborLabelCounts.indexIsUsed(neighborLabel)) {
            neighborLabelCounts.insert(neighborLabel, 0);
        }
        ++neighborLabelCounts[neighborLabel];
    });

    // Evaluate each label
    SparseVector<double> &labelWeights = labelWeightsPerThread[tid];
    assert(labelWeights.isClean());
    for (label candidateLabel : neighborLabelCounts.insertedIndexes()) {
        count localCount = neighborLabelCounts[candidateLabel];
        double weight = localCount - alpha * (globalLabelCounts[candidateLabel] - localCount);
        labelWeights.insert(candidateLabel, weight);
    }

    // Get best label
    Aux::UniformRandomSelector selector{};
    label bestLabel = none;
    double bestWeight = -std::numeric_limits<double>::max();
    for (auto neighborLabel : labelWeights.insertedIndexes()) {
        double weight = labelWeights[neighborLabel];
        if (weight > bestWeight) {
            bestWeight = weight;
            bestLabel = neighborLabel;
            selector.reset();
        } else if (weight == bestWeight) {
            // if multiple labels have the same weight, select one randomly
            if (selector.addElement())
                bestLabel = neighborLabel;
        }
    }
    assert(bestLabel != none);

    neighborLabelCounts.reset();
    labelWeights.reset();
    return bestLabel;
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
    return this->iteration;
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
