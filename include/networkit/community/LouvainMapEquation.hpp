/*
 * LouvainMapEquation.h
 *
 * Created on: 2019-01-28
 * Author: Armin Wiebigke
 *         Michael Hamann
 *         Lars Gottesbüren
 */

#include <vector>
#include <cstddef>
#include <algorithm>
#include <random>
#include <cmath>
#include <mutex>
#include <atomic>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/auxiliary/SparseVector.hpp>
#include <networkit/community/ClusteringFunctionFactory.hpp>

namespace NetworKit {

class Spinlock {
public:
    void lock() {
        while(spinner.test_and_set(std::memory_order_acquire)) {

        }
    }

    void unlock() {
        spinner.clear(std::memory_order::memory_order_release);
    }
private:
    std::atomic_flag spinner = ATOMIC_FLAG_INIT;
};



class LouvainMapEquation : public Algorithm {
private:
    friend class LouvainMapEquationFactory;
    enum class ParallelizationType : uint8_t { RelaxMap, SynchronousLocalMoving };
public:
    explicit LouvainMapEquation(const Graph &graph, bool hierarchical = false, count maxIterations = 256,
                                bool parallel = false, ParallelizationType parallelizationType = ParallelizationType::SynchronousLocalMoving,
                                double additionalCut = 0.0, double additionalVolume = 0.0);

    void run() override;

    /**
     * Returns the result of the algorithm.
     */
    const Partition& getPartition() const;

    std::string toString() const override;

private:
    struct Move {
        node movedNode;
        double volume;
        index cacheID, originCluster, targetCluster;
        double cutUpdateToOriginCluster, cutUpdateToTargetCluster;

        Move(const node n = none, double vol = 0.0, index c = none, index cc = none, index tc = none, double cuptoc = 0.0, double cupttc = 0.0) :
                movedNode(n), volume(vol), cacheID(c), originCluster(cc), targetCluster(tc), cutUpdateToOriginCluster(cuptoc), cutUpdateToTargetCluster(cupttc) { }
    };

    struct NeighborInChunk {
        node neighbor;
        index oldCluster;
        double weightToNeighbor;
        NeighborInChunk(node n = none, index oc = none, double wtn = 0.0) : neighbor(n), oldCluster(oc), weightToNeighbor(wtn) { }
    };

    static_assert(std::is_trivially_destructible<Move>::value, "LouvainMapEquation::Move struct is not trivially destructible");
    static_assert(std::is_trivially_destructible<NeighborInChunk>::value, "LouvainMapEquation::NeighborInChunk struct is not trivially destructible");

    using NeighborCache = std::vector<NeighborInChunk>;
    using NeighborCaches = std::vector<NeighborCache>;

    const bool parallel;
    ParallelizationType parallelizationType;

    const Graph& graph;
    bool hierarchical;
    count maxIterations;

    Partition partition;
    std::vector<double> clusterCut, clusterVolume;
    const double additionalCut, additionalVolume;
    double totalCut, totalVolume;

    // for RelaxMap
    std::vector< Spinlock > locks;

    // for SLM
    Partition nextPartition;
    std::vector< SparseVector<double> > ets_neighborClusterWeights;
    std::vector< std::vector<bool> > ets_isNodeInCurrentChunk;
    std::vector< NeighborCaches > ets_neighborCaches;

#ifndef NDEBUG
    long double sumPLogPwAlpha = 0;
    long double sumPLogPClusterCut = 0;
    long double sumPLogPClusterCutPlusVol = 0;
    double plogpRel(count w);
    void updatePLogPSums();
    double mapEquation();
    void checkUpdatedCutsAndVolumesAgainstRecomputation();
#endif

    count localMoving(std::vector<node>& nodes, count iteration);

    count synchronousLocalMoving(std::vector<node>& nodes, count iteration);

    template<bool parallel, bool synchronous>
    bool tryLocalMove(node u, SparseVector<double>& neighborClusterWeights, index& cacheID, std::vector<NeighborInChunk>& cachedNeighbors, std::vector<Move>& moves, std::vector<bool>& isNodeInCurrentChunk);

    template<bool parallel>
    bool performMove(node u, double degree, double loopWeight, node currentCluster, node targetCluster, double weightToTarget, double weightToCurrent);

    void aggregateAndApplyCutAndVolumeUpdates(std::vector<Move>& moves, NeighborCaches& neighborCaches);

    void calculateInitialClusterCutAndVolume();

    void runHierarchical();

    /**
    * Calculate the change in the map equation if the node is moved from its current cluster to the target cluster.
    * To simplify the calculation, we remove terms that are constant for all target clusters. As a result, "moving" the
    * node to its current cluster gives a value != 0, although the complete map equation would not change.
    */
    double fitnessChange(node, double degree, double loopWeight,
                         node currentCluster, node targetCluster,
                         double weightToTarget, double weightToCurrent, double totalCutCurrently) {
        const double cutTarget = clusterCut[targetCluster];
        const double volTarget = clusterVolume[targetCluster];
        const double cutDifferenceCurrent = 2 * weightToCurrent - degree + 2 * loopWeight;
        double totalCutNew, targetClusterCutNew, targetClusterCutCurrent, targetCutPlusVolumeNew, targetCutPlusVolumeCurrent;
        if (currentCluster != targetCluster) {
            double cutDifferenceTarget = degree - 2 * weightToTarget - 2 * loopWeight;

            totalCutNew = totalCutCurrently + cutDifferenceCurrent + cutDifferenceTarget;
            targetClusterCutNew = cutTarget + cutDifferenceTarget;
            targetClusterCutCurrent = cutTarget;
            targetCutPlusVolumeNew = cutTarget + cutDifferenceTarget + volTarget + degree;
            targetCutPlusVolumeCurrent = cutTarget + volTarget;
        } else {
            totalCutNew = totalCutCurrently;
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

        normalizeAndPLogP(totalCutNew);
        normalizeAndPLogP(targetClusterCutNew);
        normalizeAndPLogP(targetClusterCutCurrent);
        normalizeAndPLogP(targetCutPlusVolumeNew);
        normalizeAndPLogP(targetCutPlusVolumeCurrent);

        return totalCutNew + ((targetCutPlusVolumeNew - targetCutPlusVolumeCurrent) - (2 * (targetClusterCutNew - targetClusterCutCurrent)));
    }

};

class LouvainMapEquationFactory : public ClusteringFunctionFactory {
public:
    explicit LouvainMapEquationFactory(bool hierarchical = false, count maxIterations = 256, std::string parallelization = "none");

    ClusteringFunction getFunction() const override;

private:
    bool hierarchical;
    count maxIterations;
    std::string parallelization;
};

}

