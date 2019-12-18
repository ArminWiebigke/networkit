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
    std::vector<std::vector<int64_t>> globalLabelCountChangePerThread(
            parallel ? numThreads : 0, std::vector<int64_t>(nodeIdBound, 0));
    std::vector<u_int8_t> activeNodes(nodeIdBound);
    G.forNodes([&](node u) { activeNodes[u] = true; });
    std::vector<std::vector<u_int8_t>> nextActiveNodesPerThread(
            parallel ? numThreads : 0, std::vector<u_int8_t>(nodeIdBound, false));

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
            auto &globalLabelCountChange = parallel ? globalLabelCountChangePerThread[threadId]
                                                    : globalLabelCounts;
            auto &nextActiveNodesThisThread = parallel ? nextActiveNodesPerThread[threadId]
                                                       : activeNodes;
            auto evaluateNode = [&](node u) {
                if (!activeNodes[u] || G.degree(u) == 0)
                    return;
                activeNodes[u] = false;

                label bestLabel = calculateBestLabel(u);

                // Update labels
                label currentLabel = result.subsetOf(u);
                if (currentLabel != bestLabel) {
                    ++localUpdatedNodes;
                    nextPartition.moveToSubset(bestLabel, u);
                    --globalLabelCountChange[currentLabel];
                    ++globalLabelCountChange[bestLabel];
                    nextActiveNodesThisThread[u] = true;
                    G.forNeighborsOf(u, [&](node w) {
                        nextActiveNodesThisThread[w] = true;
                    });
                }
            };

            if (parallel) {
#pragma omp for nowait
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
            count chunkSize = (nodeIdBound - 1) / numThreads + 1;
            chunkSize = std::max(chunkSize, (count) 128);
            const count numberOfChunks = 1 + (nodeIdBound - 1) / chunkSize;
            assert(chunkSize == 128 || numberOfChunks == numThreads);
            std::vector<index> chunkBorders = Aux::Parallel::Chunking::getChunkBorders(
                    nodeIdBound, numberOfChunks);
#pragma omp parallel for schedule(static, 1)
            for (index chunk = 0; chunk < numberOfChunks; ++chunk) {
                for (index tid = 0; tid < numThreads; ++tid) {
                    for (index u = chunkBorders[chunk]; u < chunkBorders[chunk + 1]; ++u) {
                        // for label u
                        globalLabelCounts[u] += globalLabelCountChangePerThread[tid][u];
                        globalLabelCountChangePerThread[tid][u] = 0;
                        // for node u
                        activeNodes[u] |= nextActiveNodesPerThread[tid][u];
                        nextActiveNodesPerThread[tid][u] = false;
                    }
                }
            }
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
