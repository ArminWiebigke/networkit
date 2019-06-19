/*
 * EgoNetPartition.cpp
 *
 * Created: 2019-06-17
 * Author: Armin Wiebigke
 */

#include <memory>

#include "EgoNetPartition.h"
#include "../auxiliary/Timer.h"
#include "../coarsening/ParallelPartitionCoarsening.h"
#include "../structures/NodeMapping.h"
#include "../oslom/Stochastics.h"
#include "ExtendSignificance.h"
#include "ExtendEdges.h"

namespace NetworKit {

EgoNetPartition::EgoNetPartition(const EgoNetData &egoNetData,
                                 const PartitionFunction &partitionFunction)
		: CommunityDetectionAlgorithm(egoNetData.G),
		  directedG(egoNetData.directedG),
		  egoGraph(egoNetData.egoGraph),
		  egoMapping(egoNetData.egoMapping),
		  egoNode(egoNetData.egoNode),
		  partitionFunction(partitionFunction),
		  parameters(egoNetData.parameters),
		  groundTruth(egoNetData.groundTruth),
		  egoNetData(egoNetData) {

}

bool EgoNetPartition::isParallel() const {
	return false;
}

std::string EgoNetPartition::toString() const {
	return "EgoNetPartition";
}


void EgoNetPartition::run() {
	Aux::Timer timer;
	timer.start();

	// TODO: EgoNet Aufbau woanders hin?
	/******************************************************************************************
	 **                                Add Neighbor Nodes                                    **
	 ******************************************************************************************/
	INFO("Add neighbors");

	if (parameters.at("addEgoNode") == "Yes") {
		egoMapping.addNode(egoNode);
		egoGraph.addNode();
	}
	// Add neighbors
	G.forEdgesOf(egoNode, [&](node, node v) {
		egoMapping.addNode(v);
	});
	addTime(timer, "1    Find nodes");


	/******************************************************************************************
	 **                             Triangle Search for Edges                                **
	 ******************************************************************************************/
	INFO("Add edges");
	// Find all triangles and add the edges to the egoGraph
	G.forEdgesOf(egoNode, [&](node, node v, edgeweight weight1) {
		if (parameters.at("addEgoNode") == "Yes") {
			egoGraph.addEdge(egoMapping.local(egoNode), egoMapping.local(v), weight1);
		}
		directedG.forEdgesOf(v, [&](node, node w, edgeweight weight2) {
			if (egoMapping.isMapped(w)) {
				// we have found a triangle u-v-w
				egoGraph.addEdge(egoMapping.local(v), egoMapping.local(w), weight2);
			}
		});
	});
	addTime(timer, "2    Neighbor Triangle Search");


	/******************************************************************************************
	 **                               Extend and Partition                                   **
	 ******************************************************************************************/
	INFO("Extend EgoNet");
	Graph egoGraphBase{egoGraph};
	NodeMapping nodeMappingBase{egoMapping};
	count extendIterations = std::stoi(parameters.at("extendPartitionIterations"));
	INFO("Extend for " + std::to_string(extendIterations) + " iterations");
	std::string extendStrategy = parameters.at("extendStrategy");
	if (extendIterations < 1) {
		partitionEgoNet();
		addTime(timer, "4    Partition EgoNet");
	}
	for (count i = 0; i < extendIterations; ++i) {
		if (i > 0) {
			extendStrategy = parameters.at("extendStrategySecond");
			egoGraph = egoGraphBase;
			egoMapping = nodeMappingBase;
		}
		INFO("Extend ego-net with strategy " + extendStrategy);

		extendEgoNet(extendStrategy);
		addTime(timer, "3    Extend EgoNet");

		partitionEgoNet();
		addTime(timer, "4    Partition EgoNet");
	}

	hasRun = true;
}

void EgoNetPartition::partitionEgoNet() {
	if (parameters.at("partitionFromGroundTruth") == "Yes")
		result = createGroundTruthPartition();
	else if (egoGraph.numberOfEdges() > 0) {
		result = partitionFunction(egoGraph);
	} else {
		result = Partition(egoGraph.upperNodeIdBound());
		result.allToSingletons();
	}
	result.compact();
}

Partition
EgoNetPartition::createGroundTruthPartition() const {
	auto truthComms = groundTruth.subsetsOf(egoNode);
	Partition part{egoGraph.upperNodeIdBound()};
	egoGraph.forNodes([&](node v) {
		auto comms = groundTruth.subsetsOf(egoMapping.global(v));
		std::vector<node> overlap(truthComms.size());
		std::set_intersection(truthComms.begin(), truthComms.end(), comms.begin(), comms.end(),
		                      overlap.begin());
		if (!overlap.empty()) {
			part.addToSubset(overlap[0], v);
		}
	});
	return part;
}

void
EgoNetPartition::extendEgoNet(const std::string &extendStrategy) {
	if (extendStrategy == "none")
		return;

	Aux::Timer timer;
	timer.start();
	const count directNeighborsCnt = egoMapping.nodeCount();
	bool useBasePartition = result.numberOfSubsets() > 0;
	if (!useBasePartition) {
		result = Partition(egoGraph.upperNodeIdBound());
		result.allToOnePartition();
	}
	addTime(timer, "31    Setup");

//	std::set<node> foundNeighbors;
//	std::vector<node> directNeighbors = egoMapping.globalNodes();
//	auto isDirectNeighbor = [&](node x) {
//		return egoMapping.isMapped(x);
//	};
//	for (node v : directNeighbors) {
//		directedG.forEdgesOf(v, [&](node v, node w, edgeweight weight) {
//			if (!isDirectNeighbor(w) && w != egoNode) {
//				foundNeighbors.insert(w);
//			}
//		});
//	}


	/**********************************************************************************************
	 **                           Get node candidates with scores                                **
	 **********************************************************************************************/
	std::vector<std::pair<node, double>> nodeScores; // node and its score
	std::unique_ptr<ExtendScore> extendScore;
	if (extendStrategy == "edgeScore") {
		extendScore.reset(new ExtendEdges(egoNetData));
	} else if (extendStrategy == "significance") {
		extendScore.reset(new ExtendSignificance(egoNetData, result));
	} else
		throw std::runtime_error(extendStrategy
		                         + " is not a valid strategy to extend the Ego-Net!");
	extendScore->run();
	nodeScores = extendScore->getScores();

//	if (parameters.at("onlyDirectedCandidates") == "Yes") {
//		// Remove candidates that can not be found over the directed graph
//		auto newEnd = std::remove_if(nodeScores.begin(), nodeScores.end(),
//		                             [&](std::pair<node, double> pair) {
//			                             return foundNeighbors.count(pair.first) == 0;
//		                             });
//		nodeScores.resize(newEnd - nodeScores.begin());
//	}
	addTime(timer, "33    Get candidates");

#ifndef NDEBUG
	std::set<node> candidates;
	for (auto pair : nodeScores) {
		node v = pair.first;
		assert(candidates.count(v) == 0);
		candidates.insert(v);
	}
	assert(candidates.size() == nodeScores.size());
#endif


	/**********************************************************************************************
	 **                                  Add nodes to ego-net                                    **
	 **********************************************************************************************/
	double addNodesFactor = std::stod(parameters.at("addNodesFactor"));
	double addNodesExponent = std::stod(parameters.at("addNodesExponent"));
	count extendNodeCnt = std::ceil(
			addNodesFactor * std::pow(egoGraph.numberOfNodes(), addNodesExponent));

	// Take the nodes with the best scores
	std::sort(nodeScores.begin(), nodeScores.end(),
	          [](std::pair<node, double> a, std::pair<node, double> b) {
		          return a.second > b.second;
	          });
	if (nodeScores.size() > extendNodeCnt)
		nodeScores.resize(extendNodeCnt);
	addTime(timer, "34    Sort candidates");
	for (auto pair : nodeScores) {
		if (pair.second <= 0)
			break;
		egoMapping.addNode(pair.first);
		egoGraph.addNode();
	}
	addTime(timer, "35    Add nodes");


	/**********************************************************************************************
	 **                                  Add edges to ego-net                                    **
	 **********************************************************************************************/
	// TODO: We already looked at the edges from neighbors to neighbors of neighbors, store these
	//  for each candidate. Then we only have to add edges between neighbors of neighbors.
	const bool discardNeighEdges = parameters.at("edgesBetweenNeigNeig") != "Yes";
	for (node v : egoMapping.globalNodes()) {
		directedG.forEdgesOf(v, [&](node, node w, edgeweight weight) {
			if (egoMapping.isMapped(w)) {
				node v_loc = egoMapping.local(v);
				node w_loc = egoMapping.local(w);
				// Edges between direct neighbors are already in the Egonet
				if (v_loc < directNeighborsCnt && w_loc < directNeighborsCnt)
					return;

				// Discard edges between neighbors of neighbors
				if (discardNeighEdges && v_loc >= directNeighborsCnt
				    && w_loc >= directNeighborsCnt) {
					return;
				}
				egoGraph.addEdge(v_loc, w_loc, weight);
			}
		});
	}
	addTime(timer, "37    Add edges");

	count minDegree = stoi(parameters.at("minNodeDegree"));
	removeLowDegreeNodes(minDegree, directNeighborsCnt);
	addTime(timer, "39    Remove low degree nodes");
}

void EgoNetPartition::removeLowDegreeNodes(count minDegree, count directNeighborsCnt) const {
	if (stoi(parameters.at("minNodeDegree")) <= 0)
		return;
	count nodes_changed;
	do {
		nodes_changed = 0;
		egoGraph.forNodes([&](node v) {
			if (v >= directNeighborsCnt && egoGraph.degree(v) < minDegree) {
				++nodes_changed;
				egoGraph.removeNode(v);
//				egoGraph.restoreNode(v); // For direct neighbors?
			}
		});
	} while (nodes_changed > 0);
}

} /* namespace NetworKit */