/*
 * OLP.h
 *
 *  Created on: 07.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef OLP_H_
#define OLP_H_

#include "CommunityDetectionAlgorithm.h"
#include "../structures/Partition.h"
#include "../structures/Cover.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup community
 * As described in Ovelgoenne et al: An Ensemble Learning Strategy for Graph Clustering
 * Raghavan et al. proposed a label propagation algorithm for graph clustering.
 * This algorithm initializes every vertex of a graph with a unique label. Then, in iterative
 * sweeps over the set of vertices the vertex labels are updated. A vertex gets the label
 * that the maximum number of its neighbors have. The procedure is stopped when every vertex
 * has the label that at least half of its neighbors have.
 *
 */
class OLP : public Algorithm {

protected:

    Graph G;
    count maxLabels;
    count updateThreshold = 0;
    count maxIterations;
    count nIterations = 0; //!< number of iterations in last run
    std::vector<count> timing;    //!< running times for each iteration
    Cover result;


public:

    /**
     * Constructor to the label propagation community detection algorithm.
     *
     * @param[in]	G	input graph
     * @param[in]	theta	updateThreshold: number of nodes that have to be changed in each iteration so that a new iteration starts.
     */
    OLP(const Graph &G, count k = 3, count theta = none, count maxIterations = 20);

    /**
     * Run the label propagation clustering algorithm.
     */
    void run() override;

    /**
     * @return String representation of algorithm and parameters.
     */
    std::string toString() const override;


    /**
     * The algorithm runs until a number of nodes less than
     * the threshold is updated.
     *
     * @param th The threshold.
    */
    virtual void setUpdateThreshold(count th);

    /**
    * Get number of iterations in last run.
    *
    * @return The number of iterations.
    */
    virtual count numberOfIterations();


    /**
    * Get list of running times for each iteration.
    *
    * @return The list of running times in milliseconds
    */
    virtual std::vector<count> getTiming();

    Cover getCover();


};

} /* namespace NetworKit */
#endif /* OLP_H_ */
