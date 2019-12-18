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
#include <networkit/auxiliary/Parallel.hpp>

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
    globalLabelCounts.resize(G.upperNodeIdBound(), 1);
}

void LPPotts::run() {
    if (hasRun) {
        throw std::runtime_error("The algorithm has already run on the graph.");
    }

    runAlgorithm();

    hasRun = true;
}

void LPPotts::runAlgorithm() {
    const index nodeIdBound = G.upperNodeIdBound();
    count updatedNodesCount = G.numberOfNodes();
    iteration = 0;

    // Create data structure for parallel execution
    Partition secondPartition;
    if (parallel) {
        secondPartition = result;
    }
    Partition &nextPartition = parallel ? secondPartition : result;
    const count numThreads = Aux::getMaxNumberOfThreads();
    std::vector<SparseVector<int64_t>> globalLabelCountChangePerThread(
            parallel ? numThreads : 0, SparseVector<int64_t>(nodeIdBound, 0));

    SparseVector<u_int8_t> nodeIsActive(nodeIdBound);
    G.forNodes([&](node u) {
        if (G.degree(u) > 0) {
            nodeIsActive.insert(u, true);
        }
    });
    std::vector<SparseVector<u_int8_t>> nextActiveNodesPerThread(
            parallel ? numThreads : 0, SparseVector<u_int8_t>(nodeIdBound, false));

    // Propagate labels
    while ((updatedNodesCount > this->updateThreshold) && (iteration < maxIterations)) {
        // as long as enough labels change or maximum iterations reached
        iteration += 1;
        updatedNodesCount = 0;
        DEBUG("[BEGIN] LabelPropagation: iteration #", iteration);
#pragma omp parallel if (parallel)
        {
            count localUpdatedNodes = 0;
            index threadId = getThreadId();
            assert(threadId == omp_get_thread_num());
            auto &globalLabelCountChangeThisThread = globalLabelCountChangePerThread[threadId];
            auto &nextActiveNodesThisThread = parallel ? nextActiveNodesPerThread[threadId]
                                                       : nodeIsActive;
            auto evaluateNode = [&](node u) {
                if (!nodeIsActive[u])
                    return;
                nodeIsActive[u] = false;

                label bestLabel = calculateBestLabel(u);

                // Update labels
                label currentLabel = result.subsetOf(u);
                if (currentLabel != bestLabel) {
                    ++localUpdatedNodes;
                    nextPartition.moveToSubset(bestLabel, u);

                    if (parallel) {
                        globalLabelCountChangeThisThread.insertIfEmpty(currentLabel, 0);
                        globalLabelCountChangeThisThread.insertIfEmpty(bestLabel, 0);
                        --globalLabelCountChangeThisThread[currentLabel];
                        ++globalLabelCountChangeThisThread[bestLabel];
                    } else {
                        --globalLabelCounts[currentLabel];
                        ++globalLabelCounts[bestLabel];
                    }

                    nextActiveNodesThisThread.insertOrSet(u, true);
                    G.forNeighborsOf(u, [&](node w) {
                        nextActiveNodesThisThread.insertOrSet(w, true);
                    });
                }
            };

            if (parallel) {
                auto &activeNodes = nodeIsActive.insertedIndexes();
#pragma omp for nowait
                for (index i = 0; i < activeNodes.size(); ++i) {
                    evaluateNode(activeNodes[i]);
                }
            } else {
                G.forNodesInRandomOrder(evaluateNode);
            }

#pragma omp atomic update
            updatedNodesCount += localUpdatedNodes;
        }

        if (parallel) {
            nodeIsActive.reset();
            for (index tid = 0; tid < numThreads; ++tid) {
                for (node u : nextActiveNodesPerThread[tid].insertedIndexes()) {
                    assert(nextActiveNodesPerThread[tid][u]);
                    nodeIsActive.insertOrSet(u, true);
                }
                nextActiveNodesPerThread[tid].reset();

                for (label l : globalLabelCountChangePerThread[tid].insertedIndexes()) {
                    globalLabelCounts[l] += globalLabelCountChangePerThread[tid][l];
                }
                globalLabelCountChangePerThread[tid].reset();
            }

            result = nextPartition;
        }

//		this->timing.push_back(runtime.elapsedMilliseconds());
//		DEBUG("[DONE] LabelPropagation: iteration #", iteration, " - updated ", updatedNodesCount,
//		      " labels, time spent: ",
//		      runtime.elapsedTag());
    } // end while
}

label LPPotts::calculateBestLabel(node u) {
    index threadId = getThreadId();

    // Count the labels of the neighbors
    SparseVector<count> &neighborLabelCounts = neighborLabelCountsPerThread[threadId];
    assert(neighborLabelCounts.isClean());
    G.forNeighborsOf(u, [&](node w) {
        label neighborLabel = result.subsetOf(w);
        if (!neighborLabelCounts.indexIsUsed(neighborLabel)) {
            neighborLabelCounts.insert(neighborLabel, 0);
        }
        ++neighborLabelCounts[neighborLabel];
    });

    // Evaluate each label
    SparseVector<double> &labelWeights = labelWeightsPerThread[threadId];
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

index LPPotts::getThreadId() const {
    return parallel ? omp_get_thread_num() : 0;
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
