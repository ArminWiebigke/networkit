/*
 * ParallelPartitionCoarsening.h
 *
 *  Created on: 03.07.2014
 *      Author: cls
 */

#ifndef PARALLELPARTITIONCOARSENING_H_
#define PARALLELPARTITIONCOARSENING_H_

#include "../Globals.h"
#include "GraphCoarsening.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 * @ingroup coarsening
 */
class ParallelPartitionCoarsening: public GraphCoarsening {
public:
	ParallelPartitionCoarsening(const Graph& G, const Partition& zeta, bool useGraphBuilder = true, bool parallel = true);

	void run() override;

private:
	const Partition& zeta;
	bool useGraphBuilder;
	bool parallel;
};

} /* namespace NetworKit */

#endif /* PARALLELPARTITIONCOARSENING_H_ */
