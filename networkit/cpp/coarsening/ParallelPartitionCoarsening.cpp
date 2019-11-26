/*
 * ParallelPartitionCoarsening.cpp
 *
 *  Created on: 28.01.2014
 *      Author: cls
 */

#include <numeric>
#include <cassert>


#include "ParallelPartitionCoarsening.h"
#include <omp.h>
#include "../graph/GraphBuilder.h"
#include "../auxiliary/Timer.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

ParallelPartitionCoarsening::ParallelPartitionCoarsening(const Graph& G,
		const Partition& zeta, bool useGraphBuilder, bool parallel)
		: GraphCoarsening(G), zeta(zeta), useGraphBuilder(useGraphBuilder),
		  parallel(parallel) {

}

void ParallelPartitionCoarsening::run() {
	Partition nodeMapping = zeta;
	nodeMapping.compact(nodeMapping.upperBound() <= G.upperNodeIdBound());
	index numParts = nodeMapping.upperBound();
	
	// Leave out parallel counting sort for now as it requires some more setup.
	std::vector<index> partBegin(numParts + 2, 0);
	std::vector<node> nodesSortedByPart(G.numberOfNodes());
	G.forNodes([&](const node u) {
		partBegin[ nodeMapping[u] + 2 ]++;
	});
	std::partial_sum(partBegin.begin(), partBegin.end(), partBegin.begin());
	G.forNodes([&](const node u) {
		nodesSortedByPart[ partBegin[ nodeMapping[u] + 1 ]++ ] = u;
	});
	
	Gcoarsened = Graph(numParts, true, false);
	
	if (!parallel) {
	
		std::vector<edgeweight> incidentPartWeights(numParts, 0.0);
		std::vector<node> incidentParts;
		incidentParts.reserve(numParts);
		
		count numEdges = 0;
		count numSelfLoops = 0;
		for (node su = 0; su < numParts; ++su) {
			for (index i = partBegin[su]; i < partBegin[su+1]; ++i) {
				node u = nodesSortedByPart[i];
				G.forNeighborsOf(u, [&](node v, edgeweight ew) {
					const node sv = nodeMapping[v];
					if (sv != su || u >= v) {
						if (incidentPartWeights[sv] == 0.0) {
							incidentParts.push_back(sv);
						}
						incidentPartWeights[sv] += ew;
					}
				});
			}
			
			numEdges += incidentParts.size();
			if (incidentPartWeights[su] != 0.0) {
				numSelfLoops += 1;
				numEdges -= 1;
			}
			for (node sv : incidentParts) {
				Gcoarsened.addHalfEdge(su, sv, incidentPartWeights[sv]);
				incidentPartWeights[sv] = 0.0;
			}
			incidentParts.clear();
		}
		
		Gcoarsened.m = numEdges / 2 + numSelfLoops;
		Gcoarsened.storedNumberOfSelfLoops = numSelfLoops;

	} else {
		int t = omp_get_max_threads();
		std::vector<count> numEdges(t, 0);
		std::vector<count> numSelfLoops(t, 0);
		
		
		#pragma omp parallel
		{
			std::vector<edgeweight> incidentPartWeights(numParts, 0.0);
			std::vector<node> incidentParts;
			incidentParts.reserve(numParts);
			int tid = omp_get_thread_num();

			#pragma omp for schedule(guided)
			for (node su = 0; su < numParts; ++su) {
				for (index i = partBegin[su]; i < partBegin[su+1]; ++i) {
					node u = nodesSortedByPart[i];
					G.forNeighborsOf(u, [&](node v, edgeweight ew) {
						const node sv = nodeMapping[v];
						if (sv != su || u >= v) {
							if (incidentPartWeights[sv] == 0.0) {
								incidentParts.push_back(sv);
							}
							incidentPartWeights[sv] += ew;
						}
					});
				}
				
				numEdges[tid] += incidentParts.size();
				if (incidentPartWeights[su] != 0.0) {
					numSelfLoops[tid] += 1;
					numEdges[tid] -= 1;
				}
				for (node sv : incidentParts) {
					Gcoarsened.addHalfEdge(su, sv, incidentPartWeights[sv]);
					incidentPartWeights[sv] = 0.0;
				}
		
				incidentParts.clear();
			}
			
		}
		
		Gcoarsened.m = 0;
		Gcoarsened.storedNumberOfSelfLoops = 0;
		for (size_t tid = 0; tid < numEdges.size(); ++tid) {
			Gcoarsened.m += numEdges[tid];
			Gcoarsened.storedNumberOfSelfLoops += numSelfLoops[tid];
		}
		Gcoarsened.m =  Gcoarsened.m / 2 + Gcoarsened.storedNumberOfSelfLoops;
	}

	this->nodeMapping = nodeMapping.moveVector();
	hasRun = true;
}

} /* namespace NetworKit */
