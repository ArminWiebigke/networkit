/*
 * EgoSplitting.h
 *
 * Created: 2018-12-11
 * Author: Armin Wiebigke
 */

#include <omp.h>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "EgoSplitting.h"
#include "../structures/Partition.h"
#include "../components/ConnectedComponents.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/SignalHandling.h"
#include "../auxiliary/Timer.h"
#include "../coarsening/ParallelPartitionCoarsening.h"
#include "../graph/RandomMaximumSpanningForest.h"
#include "PLM.h"
#include "EgoNetPartition.h"
#include "../oslom/Stochastics.h"
#include "../auxiliary/ParseString.h"

#define true_or_throw(cond, msg) if (!cond) throw std::runtime_error(msg)
#define W(x) #x << "=" << x << ", "

namespace NetworKit {

EgoSplitting::EgoSplitting(const Graph &G)
		: G(G) {
	PartitionFunction clusterAlgo = [](const Graph &G) {
		PLM algo(G, false, 1.0, "none");
		algo.run();
		return algo.getPartition();
	};
	localClusterAlgo = clusterAlgo;
	globalClusterAlgo = clusterAlgo;
	init();
}

EgoSplitting::EgoSplitting(const Graph &G,
                           PartitionFunction clusterAlgo)
		: EgoSplitting(G, clusterAlgo, std::move(clusterAlgo)) {
}

EgoSplitting::EgoSplitting(const Graph &G,
                           PartitionFunction localClusterAlgo,
                           PartitionFunction globalClusterAlgo)
		: G(G),
		  localClusterAlgo(std::move(localClusterAlgo)),
		  globalClusterAlgo(std::move(globalClusterAlgo)) {
	init();
}

void EgoSplitting::init() {
	hasRun = false;
	egoNets.resize(G.upperNodeIdBound());
	personaEdges.resize(G.upperNodeIdBound());
	egoNetPartitions.resize(G.upperNodeIdBound());
	egoNetPartitionCounts.resize(G.upperNodeIdBound(), 0);
	personaOffsets.resize(G.upperNodeIdBound() + 1, 0);
	directedG = AdjacencyArray(G);

	parameters["storeEgoNet"] = "No";
	parameters["addEgoNode"] = "No";
	parameters["partitionFromGroundTruth"] = "No";
	parameters["maxEgoNetsStored"] = "2000";

	// Connect Personas
	parameters["connectPersonas"] = "Yes";
	parameters["connectPersonasStrat"] = "spanning";
	parameters["personaEdgeWeightFactor"] = "1";
	parameters["normalizePersonaCut"] = "No";
	parameters["normalizePersonaWeights"] = "unweighted";
	parameters["maxPersonaEdges"] = "1";
	parameters["iterationWeight"] = "No";

	// Parameters for weighted edges
	parameters["weightFactor"] = "1";
	parameters["weightOffset"] = "0";

	// Parameters for ego-net extension
	parameters["extendOverDirected"] = "No";
	parameters["Maximum Extend Factor"] = "1";
	parameters["addNodesExponent"] = "0.8";
	parameters["edgesBetweenNeigNeig"] = "Yes";
	parameters["minNodeDegree"] = "0";
	parameters["triangleThreshold"] = "0";
	parameters["onlyDirectedCandidates"] = "Yes";
	parameters["extendDirectedBack"] = "Yes";
	parameters["Extend and Partition Iterations"] = "1";

	// Parameters for significance extension
	parameters["Extend and Partition Iterations"] = "2";
	parameters["Sig Extend Base Clustering"] = "None";
	parameters["maxSignificance"] = "0.1";
	parameters["orderedStatPos"] = "0.1";
	parameters["sortGroups"] = "No";
	parameters["useSignInterpol"] = "Yes";
	parameters["maxGroupsConsider"] = "5";
	parameters["signMerge"] = "Yes";
	parameters["useSigMemo"] = "Yes";
	parameters["minEdgesToGroupSig"] = "1";
	parameters["sigSecondRoundStrat"] = "updateCandidates";
	parameters["secondarySigExtRounds"] = "3";
	parameters["onlyCheckSignOfMaxCandidates"] = "Yes";
	parameters["Check Candidates Factor"] = "10";
	parameters["onlyUpdatedCandidates"] = "Yes";

	// Parameters for Edgess
	parameters["Extend EgoNet Strategy"] = "Edges";
	parameters["Edges Score Strategy"] = "Edges^2 / Degree";
	parameters["Edges Iterative"] = "No";

	// Test parameters
//	parameters["storeEgoNet"] = "Yes";
//	parameters["connectPersonas"] = "No";
//	parameters["processEgoNet"] = "extend";
//	parameters["Extend EgoNet Strategy"] = "Significance";
//	parameters["Edges Score Strategy"] = "score";

}

void EgoSplitting::run() {
	timings.clear();
	Aux::SignalHandler handler;
	Aux::Timer timer;
	timer.start();
	if (Aux::stringToDouble("0.3") == 0)
		throw std::runtime_error("Can't convert numbers because of wrong locale!");
	if (G.numberOfSelfLoops() > 0)
		throw std::runtime_error("No self-loops allowed!");
//	addTime(timer, "0  Setup");

	INFO("create EgoNets");
	createEgoNets();
	addTime(timer, "1  create EgoNets");
	handler.assureRunning();

	INFO("split into Personas");
	splitIntoPersonas();
	addTime(timer, "2  split Personas");
	handler.assureRunning();

	INFO("connect Personas");
	connectPersonas();
	addTime(timer, "3  connect Personas");
	handler.assureRunning();

	INFO("create Persona Clustering");
	createPersonaClustering();
	addTime(timer, "4  Persona Clustering");
	handler.assureRunning();

	INFO("create Cover");
	createCover();
	addTime(timer, "5  create Cover");

	hasRun = true;
}

std::string EgoSplitting::toString() const {
	return "EgoSplitting";
}

void EgoSplitting::createEgoNets() {
	NodeMapping egoMapping(G); // Assign local IDs to the neighbors

	MemoizationTable<double> sigTable(-1.0, G.upperNodeIdBound());
	std::vector<double> nodeScores(G.upperNodeIdBound());
	std::vector<node> significantGroup(G.upperNodeIdBound(), none);
	std::vector<std::vector<count>> edgesToGroups(G.upperNodeIdBound());
	Aux::SignalHandler handler;
	Aux::Timer timer;

	G.forNodes([&](node u) {
		INFO("Create EgoNet for Node ", u, "/", G.upperNodeIdBound());
		handler.assureRunning();
		timer.start();
		count degree = G.degree(u);
		Graph egoGraph(degree, true);
		EgoNetData egoNetData{G, directedG, groundTruth, u, egoGraph, egoMapping, parameters,
						sigTable, nodeScores, significantGroup, edgesToGroups};


		INFO("Extend and partition EgoNet");
		EgoNetPartition extAndPartition{egoNetData, localClusterAlgo};
		extAndPartition.run();
		Partition egoPartition = extAndPartition.getPartition();
		addTimings(extAndPartition.getTimings(), "11");
		addTime(timer, "11    Extend and Partition EgoNet");


		INFO("Connect Personas of one node");
		if (parameters.at("connectPersonas") == "Yes")
			personaEdges[u] = connectEgoPartitionPersonas(egoGraph, egoPartition);
		addTime(timer, "13    Connect Personas");


		INFO("Store EgoNet");
		if (parameters.at("storeEgoNet") == "Yes")
			storeEgoNet(egoGraph, egoMapping, u);
		addTime(timer, "15    Copy EgoNet");


		INFO("Build EgoNet Partition Map");
		// Insert nodes into ego-net data structure
		for (node i : egoGraph.nodes())
			egoNetPartitions[u].emplace(egoMapping.global(i), egoPartition.subsetOf(i));
		assert(egoNetPartitions[u].size() == egoGraph.numberOfNodes());
		egoNetPartitionCounts[u] = egoPartition.numberOfSubsets();
		addTime(timer, "17    Store EgoNet Partition");


		INFO("Reset egoMapping");
		egoMapping.reset();
		addTime(timer, "1x    Clean up");
	});
}


void EgoSplitting::storeEgoNet(const Graph &egoGraph, const NodeMapping &egoMapping, node egoNode) {
	// Only store maxEgoNets ego-nets (expected)
	count maxEgoNets = std::stoi(parameters.at("maxEgoNetsStored"));
	double storePercnt = maxEgoNets * 1.0 / G.numberOfNodes();
	if (Aux::Random::real() > storePercnt)
		return;

	// Get EgoNet with gloabl node ids
	Graph egoNetGraph(G.upperNodeIdBound(), egoGraph.isWeighted());
	egoGraph.forEdges([&](node v, node w, edgeweight weight) {
		egoNetGraph.addEdge(egoMapping.global(v), egoMapping.global(w), weight);
	});
	for (node v : egoNetGraph.nodes()) {
		if (!egoGraph.hasNode(egoMapping.local(v)))
			egoNetGraph.removeNode(v);
	}
//	if (egoGraph.isWeighted() && egoNetGraph.numberOfNodes() >= 2)
//		egoNetGraph.addEdge(egoNetGraph.nodes()[0], egoNetGraph.nodes()[1],
//		                    0.000001);
	egoNetGraph.shrinkToFit();
	egoNets[egoNode] = egoNetGraph;
}

std::vector<EgoSplitting::Edge>
EgoSplitting::connectEgoPartitionPersonas(const Graph &egoGraph,
                                          const Partition &egoPartition) const {
	std::vector<EgoSplitting::Edge> edges;

	// Contract graph
	ParallelPartitionCoarsening coarsening{egoGraph, egoPartition};
	coarsening.run();
	Graph coarseGraph = coarsening.getCoarseGraph();
	auto nodeMapping = coarsening.getCoarseToFineNodeMapping();
	auto getPersonaIndex = [&](node u) {
		return egoPartition.subsetOf(nodeMapping[u][0]);
	};
	auto addPersonaEdge = [&](node u, node v, edgeweight w = 1) {
		index p_u = getPersonaIndex(u);
		index p_v = getPersonaIndex(v);
		if (p_u != p_v)
			edges.emplace_back(p_u, p_v, w);
	};

	// Get spanning forest size (for normalization)
	double spanSize;
	{
		RandomMaximumSpanningForest span{coarseGraph};
		span.run();
		auto spanningForest = span.getMSF();
		spanSize = spanningForest.numberOfEdges();
	}

	// Normalize edgeweights of the coarse graph
	if (parameters.at("normalizePersonaCut") == "volume") {
		coarseGraph.forEdges([&](node u, node v, edgeweight w) {
			double volume = w + coarseGraph.weight(u, u) + coarseGraph.weight(v, v);
			edgeweight newWeight = w / volume;
			coarseGraph.setWeight(u, v, newWeight);
		});
	} else if (parameters.at("normalizePersonaCut") == "density") {
		coarseGraph.forEdges([&](node u, node v, edgeweight w) {
			double possibleEdges = nodeMapping[u].size() * nodeMapping[v].size();
			double newWeight = w / possibleEdges;
			coarseGraph.setWeight(u, v, newWeight);
		});
	}

	// Insert edges between the personas
	std::string strategy = parameters.at("connectPersonasStrat");
	if (strategy == "spanning") {
		// TODO: Better weights
		// Every persona gets at most 'iterations' edges
		count iterations = std::stoi(parameters.at("maxPersonaEdges"));
		for (int i = 0; i < iterations; ++i) {
			RandomMaximumSpanningForest span{coarseGraph};
			span.run();
			auto spanningForest = span.getMSF();
			spanningForest.forEdges([&](node u, node v, edgeweight w) {
				if (parameters.at("iterationWeight") == "Yes")
					w = 1.0 / (i + 1);
				addPersonaEdge(u, v, w);
				coarseGraph.removeEdge(u, v);
			});
		}
	} else if (strategy == "maxEdge") {
		count iterations = std::stoi(parameters.at("maxPersonaEdges"));
		for (int i = 0; i < iterations; ++i) {
			std::vector<std::set<node>> edgeSets{coarseGraph.upperNodeIdBound()};
			coarseGraph.forNodes([&](node u) {
				edgeweight maxWeight = 0.0;
				node maxNode = none;
				coarseGraph.forEdgesOf(u, [&](node, node v, edgeweight w) {
					if (w > maxWeight) {
						maxWeight = v;
						maxNode = v;
					}
				});
				if (maxNode == none)
					return;
				if (u > maxNode)
					std::swap(u, maxNode);
				edgeSets[u].insert(maxNode);
			});
			for (node u = 0; u < edgeSets.size(); ++u) {
				for (auto v : edgeSets[u]) {
					addPersonaEdge(u, v, coarseGraph.weight(u, v));
					coarseGraph.removeEdge(u, v);
				}
			}
		}
	} else if (strategy == "all") {
		coarseGraph.forEdges([&](node u, node v, edgeweight w) {
			addPersonaEdge(u, v, w);
		});
	} else {
		throw std::runtime_error(strategy + " is not a valid strategy to connect personas!");
	}

	// Normalize the weights of the edges between the personas
	if (parameters.at("normalizePersonaWeights") == "spanSize") {
		edgeweight weightSum = 0.0;
		for (auto &edge : edges)
			weightSum += edge.weight;
		weightSum /= spanSize; // Sum of edge weights == spanning tree size
		for (auto &edge : edges) {
			edge.weight /= weightSum;
		}
	} else if (parameters.at("normalizePersonaWeights") == "unweighted") {
		for (auto &edge : edges)
			edge.weight = 1;
	} else if (parameters.at("normalizePersonaWeights") == "max1") {
		edgeweight maxWeight = 0.0;
		for (auto &edge : edges)
			maxWeight = std::max(maxWeight, edge.weight);
		for (auto &edge : edges)
			edge.weight /= maxWeight;
	} else if (parameters.at("normalizePersonaWeights") == "sameWeights") {
		for (auto &edge : edges)
			edge.weight = spanSize / edges.size();
	}
	return edges;
}

void EgoSplitting::splitIntoPersonas() {
	count sum = 0;
	for (index i = 0; i < G.upperNodeIdBound(); ++i) {
		personaOffsets[i] = sum;
		count numPersonas = std::max(egoNetPartitionCounts[i], (count) 1);
		sum += numPersonas;
	}
	personaOffsets[G.upperNodeIdBound()] = sum;
	personaGraph = Graph(sum, true);

}

void EgoSplitting::connectPersonas() {
	auto getPersona = [&](node u, index i) {
		return personaOffsets[u] + i;
	};

	// Connect personas of each node
	double weightFactor = Aux::stringToDouble(parameters.at("personaEdgeWeightFactor"));
	G.forNodes([&](node u) {
		for (auto edge : personaEdges[u]) {
			personaGraph.addEdge(getPersona(u, edge.firstNode), getPersona(u, edge.secondNode),
			                     edge.weight * weightFactor);
		}
	});

	// Connect personas of different nodes
	G.forEdges([&](node u, node v, edgeweight weight) {
		auto idx_u = egoNetPartitions[u].find(v);
		auto idx_v = egoNetPartitions[v].find(u);
		assert(idx_u != egoNetPartitions[u].end() && idx_v != egoNetPartitions[v].end());
		personaGraph.addEdge(getPersona(u, idx_u->second), getPersona(v, idx_v->second),
		                     weight);
	});

//	egoNetPartitions.clear();

#ifndef NDEBUG
	count internalPersonaEdges = 0;
	for (const auto &edges : personaEdges)
		internalPersonaEdges += edges.size();
	assert(personaGraph.numberOfEdges() == G.numberOfEdges() + internalPersonaEdges);
	// check that no isolated nodes were added
	auto numIsolatedNodes = [](const Graph &graph) {
		ConnectedComponents compsAlgo(graph);
		compsAlgo.run();
		auto comps = compsAlgo.getComponentSizes();
		count isolated = 0;
		for (const auto &x  : comps) {
			count num_nodes = x.second;
			if (num_nodes == 1)
				++isolated;
		}
		return isolated;
	};
	auto iso1 = numIsolatedNodes(G);
	auto iso2 = numIsolatedNodes(personaGraph);
	assert(iso1 == iso2);
#endif
}

void EgoSplitting::createPersonaClustering() {
	personaPartition = globalClusterAlgo(personaGraph);
	assert(personaPartition.upperBound() <= personaGraph.upperNodeIdBound());
}

void EgoSplitting::createCover() {
	// Create cover from persona partition
	cover = Cover{G.upperNodeIdBound()};
	personaPartition.compact();
	cover.setUpperBound(personaPartition.upperBound());
	G.forNodes([&](node u) {
		for (index i = personaOffsets[u]; i < personaOffsets[u + 1]; ++i) {
			cover.addToSubset(personaPartition.subsetOf(i), u);
		}
	});

	// Discard communities of size 4 or less
	count min_size = 5;
	std::vector<std::vector<node>> communities{cover.upperBound()};
	G.forNodes([&](node u) {
		for (index c : cover.subsetsOf(u)) {
			if (communities[c].size() < min_size)
				communities[c].push_back(u);
		}
	});
	for (index c = 0; c < communities.size(); ++c) {
		if (communities[c].size() < min_size) {
			for (node u : communities[c])
				cover.removeFromSubset(c, u);
		}
	}
}

Cover EgoSplitting::getCover() {
	return cover;
}

std::unordered_map<std::string, double> EgoSplitting::getExecutionInfo() {
	return executionInfo;
}

std::vector<Graph> EgoSplitting::getEgoNets() {
	return egoNets;
}

std::vector<std::unordered_map<node, index>> EgoSplitting::getEgoNetPartitions() {
	return egoNetPartitions;
}

void
EgoSplitting::setParameters(std::map<std::string, std::string> const &new_parameters) {
	for (auto &x : new_parameters) {
		this->parameters[x.first] = x.second;
	}
}

void EgoSplitting::setGroundTruth(const Cover &gt) {
	this->groundTruth = gt;

}

} /* namespace NetworKit */
