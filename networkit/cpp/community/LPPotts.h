/*
 * LPPotts.h
 *
 * Created on: 2019-01-14
 * Author: Armin Wiebigke
 */

#ifndef LPPOTTS_H_
#define LPPOTTS_H_

#include "CommunityDetectionAlgorithm.h"
#include "../structures/Partition.h"
#include "ClusteringFunctionFactory.h"

namespace NetworKit {

/**
 * @ingroup community
 * Label propagation algorithm using the Absolute Potts Model technique.
 *
 */
class LPPotts : public CommunityDetectionAlgorithm {

protected:

	double alpha;
	count updateThreshold = 0;
	count maxIterations;
	count nIterations = 0; //!< number of iterations in last run
	std::vector<count> timing;    //!< running times for each iteration
	bool parallelPropagation;


public:

	/**
	 * Constructor to the label propagation community detection algorithm.
	 *
	 * @param[in]	G	input graph
	 * @param[in]	theta	updateThreshold: number of nodes that have to be changed in each iteration so that a new iteration starts.
	 */
	explicit LPPotts(const Graph &G, double alpha = 0.3, count theta = none,
					 count maxIterations = 20, bool parallelPropagation = false);

	/**
	 * Constructor to the label propagation community detection algorithm.
	 *
	 * @param[in]	G	input graph
	 * @param[in]	baseClustering optional; the algorithm will start from the given clustering.
	 * @param[in]	theta	updateThreshold: number of nodes that have to be changed in each iteration so that a new iteration starts.
	 */
	LPPotts(const Graph& G, const Partition &baseClustering, double alpha = 0.3, count theta = none,
	        count maxIterations = 20, bool parallelPropagation = false);

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
};

class LPPottsFactory : public ClusteringFunctionFactory {
public:
	explicit LPPottsFactory(double alpha = 0.3, count theta = none,
	                        count maxIterations = 20, bool parallelPropagation = false);

	ClusteringFunction getFunction() const override;

private:
	double alpha;
	count theta;
	count maxIterations;
	bool parallelPropagation;
};

} /* namespace NetworKit */
#endif /* LPPOTTS_H_ */
