/*
 * EgoNetExtensionAndPartition.cpp
 *
 * Created: 2019-06-17
 * Author: Armin Wiebigke
 */

#include <memory>
#include <utility>

#include "EgoNetExtensionAndPartition.h"
#include "../../coarsening/ParallelPartitionCoarsening.h"
#include "../../structures/NodeMapping.h"
#include "../../oslom/Stochastics.h"
#include "ExtendSignificance.h"
#include "ExtendByScore.h"
#include "../../auxiliary/ParseString.h"

namespace NetworKit {

EgoNetExtensionAndPartition::EgoNetExtensionAndPartition(EgoNetData &egoNetData, node egoNode, Graph egoGraph,
                                                         PartitionFunction partitionFunction)
		: CommunityDetectionAlgorithm(egoNetData.G),
		  directedG(egoNetData.directedG),
		  egoGraph(std::move(egoGraph)),
		  egoMapping(egoNetData.egoMapping),
		  egoNode(egoNode),
		  partitionFunction(std::move(partitionFunction)),
		  parameters(egoNetData.parameters),
		  groundTruth(egoNetData.groundTruth),
		  egoNetData(egoNetData) {
}

void EgoNetExtensionAndPartition::run() {
	extendAndPartition();

	hasRun = true;
}

void EgoNetExtensionAndPartition::extendAndPartition() {
	Aux::Timer timer;
	timer.start();
	count extendIterations = extendIterationsCount();
	if (extendIterations == 0) {
		partitionEgoNet();
		addTime(timer, "4    Partition EgoNet");
		return;
	}

	auto extendAndPartitionFunc = [&](std::string const &strategy) {
		extendEgoNet(strategy);
		addTime(timer, "3    Extend EgoNet");

		partitionEgoNet();
		addTime(timer, "4    Partition EgoNet");
	};
	std::string extendStrategy = parameters.at("Extend EgoNet Strategy");
	std::string firstExtendStrategy = extendStrategy;
	if (extendStrategy == "Significance")
		firstExtendStrategy = parameters.at("Significance Base Extend");

	if (extendIterations == 1) {
		extendAndPartitionFunc(firstExtendStrategy);
	} else {
		Graph egoGraphBase = Graph(egoGraph);
		NodeMapping nodeMappingBase = NodeMapping(egoMapping);
		extendAndPartitionFunc(firstExtendStrategy);
		for (count i = 1; i < extendIterations; ++i) {
			egoGraph = egoGraphBase;
			egoMapping = nodeMappingBase;
			addTime(timer, "2    Copy EgoNet/Mapping");

			extendAndPartitionFunc(extendStrategy);
		}
	}
}

count
EgoNetExtensionAndPartition::extendIterationsCount() const {
	count extendIterations = std::stoi(parameters.at("Extend and Partition Iterations"));
	if (parameters.at("Extend EgoNet Strategy") == "Significance")
		++extendIterations; // Significance needs a base partition
	return extendIterations;
}

void EgoNetExtensionAndPartition::partitionEgoNet() {
	if (parameters.at("partitionFromGroundTruth") == "Yes")
		result = createGroundTruthPartition();
	else if (egoGraph.numberOfEdges() == 0) {
		result = Partition(egoGraph.upperNodeIdBound());
		result.allToSingletons();
	} else {
		result = partitionFunction(egoGraph);
	}
//	result.compact();
}

Partition
EgoNetExtensionAndPartition::createGroundTruthPartition() const {
	auto truthComms = groundTruth.subsetsOf(egoNode);
	Partition part{egoGraph.upperNodeIdBound()};
	egoGraph.forNodes([&](node v) {
		auto comms = groundTruth.subsetsOf(egoMapping.toGlobal(v));
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
EgoNetExtensionAndPartition::extendEgoNet(const std::string &extendStrategy) {
	if (extendStrategy == "None")
		return;
	Aux::Timer timer;
	timer.start();
	assert(egoMapping.nodeCount() == egoGraph.upperNodeIdBound());
	count directNeighborsBound = egoGraph.upperNodeIdBound();

	// Get node candidates with scores
	double addNodesFactor = Aux::stringToDouble(parameters.at("Maximum Extend Factor"));
	double addNodesExponent = Aux::stringToDouble(parameters.at("addNodesExponent"));
	count extendNodeCnt = std::ceil(
			addNodesFactor * std::pow(egoGraph.numberOfNodes(), addNodesExponent));
	std::vector<node> extendNodes; // nodes and their scores

	auto getExtendNodes = [&](ExtendEgoNetStrategy &&extendEgoNetStrategy) {
		extendEgoNetStrategy.run();
		extendNodes = extendEgoNetStrategy.getNodes();
		addTimings(extendEgoNetStrategy.getTimings(), "33");
	};
	if (extendStrategy == "Edges") { ;
		getExtendNodes(ExtendByScore(egoNetData, extendNodeCnt, egoGraph, egoNode));
	} else if (extendStrategy == "Significance") {
		assert(result.numberOfElements() >= egoGraph.numberOfNodes());
		getExtendNodes(ExtendSignificance(egoNetData, result, extendNodeCnt, egoGraph,
		                                  egoNode));
	} else {
		throw std::runtime_error(extendStrategy
		                         + " is not a valid strategy to extend the Ego-Net!");
	}
	addTime(timer, "33    Get candidates");

#ifndef NDEBUG
	std::set<node> candidates;
	for (node v: extendNodes) {
		assert(candidates.count(v) == 0);
		candidates.insert(v);
	}
#endif

	// Add nodes to ego-net
	assert(extendNodes.size() <= extendNodeCnt);
	for (node v : extendNodes) {
		egoMapping.addNode(v);
		egoGraph.addNode();
	}

	// Add edges to ego-net
	for (node v : egoMapping.globalNodes()) {
		directedG.forEdgesOf(v, [&](node, node w, edgeweight weight) {
			if (egoMapping.isMapped(w)) {
				node vLocal = egoMapping.toLocal(v);
				node wLocal = egoMapping.toLocal(w);
				// Edges between direct neighbors are already in the Egonet
				if (vLocal < directNeighborsBound && wLocal < directNeighborsBound)
					return;
				egoGraph.addEdge(vLocal, wLocal, weight);
			}
		});
	}

	count minDegree = std::stoi(parameters.at("minNodeDegree"));
	removeLowDegreeNodes(minDegree, directNeighborsBound);
	addTime(timer, "3a    Add candidates to ego-net");
}

void EgoNetExtensionAndPartition::removeLowDegreeNodes(count minDegree, count directNeighborsBound) {
	if (minDegree <= 0)
		return;
	count nodes_changed;
	do {
		nodes_changed = 0;
		egoGraph.forNodes([&](node v) {
			if (v >= directNeighborsBound && egoGraph.degree(v) < minDegree) {
				++nodes_changed;
				egoGraph.removeNode(v);
			}
		});
	} while (nodes_changed > 0);
}

bool EgoNetExtensionAndPartition::isParallel() const {
	return false;
}

std::string EgoNetExtensionAndPartition::toString() const {
	return "EgoNetExtensionAndPartition";
}

Graph EgoNetExtensionAndPartition::getExtendedEgoGraph() const {
	return egoGraph;
}

} /* namespace NetworKit */