/*
 * egosplitting.h
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
#include "../../structures/Partition.h"
#include "../../components/ConnectedComponents.h"
#include "../../auxiliary/SignalHandling.h"
#include "../../coarsening/ParallelPartitionCoarsening.h"
#include "../../graph/RandomMaximumSpanningForest.h"
#include "../PLM.h"
#include "EgoNetExtensionAndPartition.h"
#include "../LouvainMapEquation.h"
#include "../cleanup/SignificanceCommunityCleanUp.h"

namespace NetworKit {

EgoSplitting::EgoSplitting(const Graph &G, bool parallelEgoNetEvaluation)
		: EgoSplitting(G, parallelEgoNetEvaluation,
		               PLMFactory(true, 1.0, "none randomized").getFunction(),
		               LouvainMapEquationFactory(true).getFunction()
) {
	std::cout << "Default EgoSplitting" << std::endl;
}

EgoSplitting::EgoSplitting(const Graph &G, bool parallelEgoNetEvaluation, ClusteringFunction clusterAlgo)
		: EgoSplitting(G, parallelEgoNetEvaluation, clusterAlgo, std::move(clusterAlgo)) {
}

EgoSplitting::EgoSplitting(const Graph &G, bool parallelEgoNetEvaluation, ClusteringFunction localClusterAlgo,
                           ClusteringFunction globalClusterAlgo)
		: Algorithm(), G(G),
		  parallelEgoNetEvaluation(parallelEgoNetEvaluation),
		  localClusteringAlgo(std::move(localClusterAlgo)),
		  globalClusteringAlgo(std::move(globalClusterAlgo)) {
	init();
}

void EgoSplitting::init() {
	personaEdges.resize(G.upperNodeIdBound());
	egoNetPartitions.resize(G.upperNodeIdBound());
	egoNetExtendedPartitions.resize(G.upperNodeIdBound());
	egoNetPartitionCounts.resize(G.upperNodeIdBound(), 0);
	personaOffsets.resize(G.upperNodeIdBound() + 1, 0);
	directedG = LowToHighDirectedGraph(G);

	parameters["storeEgoNet"] = "No";
	parameters["partitionFromGroundTruth"] = "No";
	parameters["numEgoNetsStored"] = "2000";
	parameters["Cleanup"] = "Yes";
	parameters["maxEgoNetsPartitioned"] = "-1";

	// Connect Personas
	parameters["connectPersonas"] = "Yes";
	parameters["normalizePersonaCut"] = "No";
	parameters["connectPersonasStrat"] = "spanning";
	parameters["normalizePersonaWeights"] = "unweighted";
	parameters["iterationWeight"] = "No";

	// Parameters for ego-net extension
	parameters["Maximum Extend Factor"] = "5";
	parameters["addNodesExponent"] = "0.5";
	parameters["minNodeDegree"] = "2";
	parameters["Extend and Partition Iterations"] = "1";
	parameters["Extend EgoNet Strategy"] = "Edges";

	// Parameters for Edges
	parameters["Edges Score Strategy"] = "Edges pow 2 div Degree";

	// Parameters for significance extension
//	parameters["Extend and Partition Iterations"] = "2";
	parameters["Significance Base Extend"] = "None";
	parameters["maxSignificance"] = "0.1";
	parameters["sortGroups"] = "Significance";
	parameters["maxGroupsConsider"] = "99";
	parameters["signMerge"] = "Yes";
	parameters["useSigMemo"] = "No";
	parameters["minEdgesToGroupSig"] = "1";
	parameters["secondarySigExtRounds"] = "99";
	parameters["onlyCheckSignOfMaxCandidates"] = "Yes";
	parameters["Check Candidates Factor"] = "10";
	parameters["onlyUpdatedCandidates"] = "Yes";
}

void EgoSplitting::run() {
	if (hasRun)
		throw std::runtime_error("Algorithm has already been run!");
	if (G.numberOfSelfLoops() > 0)
		throw std::runtime_error("No self-loops allowed!");
	assert(timingsEmpty());
	Aux::SignalHandler handler;
	Aux::Timer timer;
	timer.start();

	INFO("create EgoNets");
	createEgoNets();
	addTime(timer, "1  Create EgoNets");
	handler.assureRunning();

	INFO("split into Personas");
	splitIntoPersonas();
	addTime(timer, "2  Split Personas");
	handler.assureRunning();

	INFO("connect Personas");
	connectPersonas();
	addTime(timer, "3  Connect Personas");
	handler.assureRunning();

	INFO("create Persona Clustering");
	createPersonaClustering();
	addTime(timer, "4  Persona Clustering");
	handler.assureRunning();

	INFO("create Cover");
	createCover();
	addTime(timer, "5  Create Cover");

	INFO("clean up Cover");
	cleanUpCover();
	addTime(timer, "6  Cleanup Cover");

	hasRun = true;
}


void EgoSplitting::createEgoNets() {
	int maxEgoNetsPartitioned = std::stoi(parameters.at("maxEgoNetsPartitioned"));
#pragma omp parallel if (parallelEgoNetEvaluation)
	{
		Aux::SignalHandler signalHandler;
		Aux::Timer timer, totalTimer;
		timer.start();
		totalTimer.start();
		NodeMapping egoMapping(G); // Assign local IDs to the neighbors
		MemoizationTable<double> sigTable(-1.0, 0);
		SparseVector<double> nodeScores(0);
		SparseVector<node> significantGroup(0, none);
		SparseVector<std::vector<count>> edgesToGroups(0);
		StochasticSignificance stochasticSignificance(0);
		EgoNetData egoNetData{G, directedG, groundTruth, egoMapping, parameters, sigTable, nodeScores,
		                      significantGroup, edgesToGroups, stochasticSignificance};
		addTime(timer, "11    Data Setup");

#pragma omp for
		for (omp_index egoNode = 0; egoNode < static_cast<omp_index>(G.upperNodeIdBound()); ++egoNode) {
			if (!signalHandler.isRunning()) continue;

			if (!G.hasNode(egoNode) || G.degree(egoNode) == 0)
				continue;

			if (G.degree(egoNode) == 1 || (maxEgoNetsPartitioned != -1 && egoNode > maxEgoNetsPartitioned)) {
				egoNetPartitionCounts[egoNode] = 1;
				G.forNeighborsOf(egoNode, [&](node neighbor) {
					egoNetPartitions[egoNode].emplace(neighbor, 0);
				});
				continue;
			}

//		G.parallelForNodes([&](node egoNode) {
			INFO("Create EgoNet for Node ", egoNode, "/", G.upperNodeIdBound());

			Graph egoGraph(G.degree(egoNode), true);
			// Find neighbors == nodes of the ego-net
			G.forEdgesOf(egoNode, [&](node, node v) {
				egoMapping.addNode(v);
			});
			// Find all triangles and add the edges to the egoGraph
			G.forEdgesOf(egoNode, [&](node, node v, edgeweight weight1) {
				directedG.forEdgesOf(v, [&](node, node w, edgeweight weight2) {
					if (egoMapping.isMapped(w)) {
						// we have found a triangle u-v-w
						egoGraph.addEdge(egoMapping.toLocal(v), egoMapping.toLocal(w), weight2);
					}
				});
			});
			addTime(timer, "13    Build EgoNet");

			DEBUG("Extend and partition");
			EgoNetExtensionAndPartition extAndPartition(egoNetData, egoNode, egoGraph,
			                                            localClusteringAlgo);
			extAndPartition.run();
			const Partition& egoPartition = extAndPartition.getPartition();
			addTimings(extAndPartition.getTimings(), "15");
			addTime(timer, "15    Extend and Partition EgoNet");

			if (parameters.at("storeEgoNet") == "Yes") { // only for analysis
				Graph extendedEgoGraph = extAndPartition.getExtendedEgoGraph();
				storeEgoNet(extendedEgoGraph, egoMapping, egoNode);
				// Store ego-net partition with extended nodes
				extendedEgoGraph.forNodes([&](node i) {
					egoNetExtendedPartitions[egoNode].emplace(egoMapping.toGlobal(i),
					                                          egoPartition.subsetOf(i));
				});
			}
			addTime(timer, "18    Store EgoNet");

			DEBUG("Store ego-net");
//			assert(egoNetExtendedPartitions[egoNode].size() == extendedEgoGraph.numberOfNodes());
			// Remove nodes that are not directed neighbors (they were added by the ego-net extension)
			Partition directNeighborPartition(egoGraph.upperNodeIdBound());
			directNeighborPartition.setUpperBound(egoPartition.upperBound());
			egoGraph.forNodes([&](node v) {
				directNeighborPartition.addToSubset(egoPartition.subsetOf(v), v);
			});
			directNeighborPartition.compact(true);
			egoNetPartitionCounts[egoNode] = directNeighborPartition.numberOfSubsets();
//			assert(egoNetPartitionCounts[egoNode] == directNeighborPartition.upperBound());
			G.forNeighborsOf(egoNode, [&](node i) {
				egoNetPartitions[egoNode].emplace(i, directNeighborPartition.subsetOf(
						egoMapping.toLocal(i)));
			});
			addTime(timer, "19    Store EgoNet Partition");

			if (parameters.at("connectPersonas") == "Yes")
				personaEdges[egoNode] = connectEgoPartitionPersonas(egoGraph, directNeighborPartition);
			addTime(timer, "16    Connect Ego-Personas");

			egoMapping.reset();
			addTime(timer, "1x    Clean up");
		}
		addTime(totalTimer, "10    EgoNet Sum");
	}
}

void EgoSplitting::storeEgoNet(const Graph &egoGraph, const NodeMapping &egoMapping, node egoNode) {
	// Only store a given maximum of ego-nets (expected)
	count maxEgoNets = std::stoi(parameters.at("numEgoNetsStored"));
	double storeChance = (double) maxEgoNets / G.numberOfNodes();
	if (Aux::Random::real() > storeChance)
		return;

	DEBUG("Store ego-net of node ", egoNode);
	// Get EgoNet with global node ids
	std::vector<WeightedEdge> edges;
	edges.reserve(egoGraph.numberOfEdges());
	// Add loops to retain isolated nodes
	egoGraph.forNodes([&](node u) {
		node globalId = egoMapping.toGlobal(u);
		edges.emplace_back(globalId, globalId, defaultEdgeWeight);
	});
	egoGraph.forEdges([&](node u, node v, edgeweight weight) {
		edges.emplace_back(egoMapping.toGlobal(u), egoMapping.toGlobal(v), weight);
	});
	egoNets[egoNode] = edges;
}

std::vector<WeightedEdge>
EgoSplitting::connectEgoPartitionPersonas(const Graph &egoGraph,
                                          const Partition &egoPartition) const {
	INFO("Connect personas of one node");
	std::vector<WeightedEdge> edges;

	// Contract graph
	ParallelPartitionCoarsening coarsening{egoGraph, egoPartition, true, false};
	coarsening.run();
	const Graph& coarseGraph = coarsening.getCoarseGraph();
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
	const std::string normalizePersonaWeights = parameters.at("normalizePersonaWeights");
	if (normalizePersonaWeights == "spanSize" || normalizePersonaWeights == "sameWeights") {
		RandomMaximumSpanningForest span{coarseGraph};
		span.run();
		auto spanningForest = span.getMSF();
		spanSize = spanningForest.numberOfEdges();
	}

	// Normalize edgeweights of the coarse graph
	std::string normalizePersonaCut = parameters.at("normalizePersonaCut");
	if (normalizePersonaCut == "volume") {
		throw std::runtime_error("Currently unsupported");
		/*
		coarseGraph.forEdges([&](node u, node v, edgeweight w) {
			double volume = w + coarseGraph.weight(u, u) + coarseGraph.weight(v, v);
			edgeweight newWeight = w / volume;
			coarseGraph.setWeight(u, v, newWeight);
		});
		*/
	} else if (normalizePersonaCut == "density") {
		throw std::runtime_error("Currently unsupported");
		/*
		coarseGraph.forEdges([&](node u, node v, edgeweight w) {
			double possibleEdges = (double) nodeMapping[u].size() * nodeMapping[v].size();
			double newWeight = w / possibleEdges;
			coarseGraph.setWeight(u, v, newWeight);
		});
		*/
	}

	// Insert edges between the personas
	const std::string strategy = parameters.at("connectPersonasStrat");
	if (strategy == "spanning") {
		RandomMaximumSpanningForest span{coarseGraph};
		span.run();
		const Graph& spanningForest = span.getMSF();
		spanningForest.forEdges([&](node u, node v, edgeweight w) {
			addPersonaEdge(u, v, w);
		});
	} else if (strategy == "all") {
		coarseGraph.forEdges([&](node u, node v, edgeweight w) {
			addPersonaEdge(u, v, w);
		});
	} else {
		throw std::runtime_error(strategy + " is not a valid strategy to connect personas!");
	}

	// Normalize the weights of the edges between the personas
	if (normalizePersonaWeights == "spanSize") {
		edgeweight weightSum = 0.0;
		for (auto &edge : edges)
			weightSum += edge.weight;
		weightSum /= spanSize; // Sum of edge weights == spanning tree size
		for (auto &edge : edges) {
			edge.weight /= weightSum;
		}
	} else if (normalizePersonaWeights == "unweighted") {
		for (auto &edge : edges)
			edge.weight = 1;
	} else if (normalizePersonaWeights == "max1") {
		edgeweight maxWeight = 0.0;
		for (auto &edge : edges)
			maxWeight = std::max(maxWeight, edge.weight);
		for (auto &edge : edges)
			edge.weight /= maxWeight;
	} else if (normalizePersonaWeights == "sameWeights") {
		for (auto &edge : edges)
			edge.weight = spanSize / edges.size();
	}
	return edges;
}

void EgoSplitting::splitIntoPersonas() {
	count sum = 0;
	G.forNodes([&](node u) {
		personaOffsets[u] = sum;
		count numPersonas = std::max(egoNetPartitionCounts[u], (count) 1);
		sum += numPersonas;
	});
	personaOffsets[G.upperNodeIdBound()] = sum;
	personaGraph = Graph(sum, (parameters.at("normalizePersonaWeights") != "unweighted" || G.isWeighted()));
}

void EgoSplitting::connectPersonas() {
	auto getPersona = [&](node u, index i) {
		assert(i < egoNetPartitionCounts[u]);
		return personaOffsets[u] + i;
	};

	// Connect personas of each node
	G.forNodes([&](node u) {
		for (auto edge : personaEdges[u]) {
			personaGraph.addEdge(getPersona(u, edge.u), getPersona(u, edge.v),
			                     edge.weight);
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
	personaPartition = globalClusteringAlgo(personaGraph);
}

void EgoSplitting::createCover() {
	// Create cover from persona partition
	resultCover = Cover(G.upperNodeIdBound());
	personaPartition.compact();
	resultCover.setUpperBound(personaPartition.upperBound());
	G.forNodes([&](node u) {
		for (index i = personaOffsets[u]; i < personaOffsets[u + 1]; ++i) {
			resultCover.addToSubset(personaPartition.subsetOf(i), u);
		}
	});
}

void EgoSplitting::discardSmallCommunities() {// Discard communities of size 4 or less
	count min_size = 5;
	std::vector<std::vector<node>> communities{resultCover.upperBound()};
	G.forNodes([&](node u) {
		for (index c : resultCover.subsetsOf(u)) {
			if (communities[c].size() < min_size)
				communities[c].push_back(u);
		}
	});
	for (index c = 0; c < communities.size(); ++c) {
		if (communities[c].size() < min_size) {
			for (node u : communities[c])
				resultCover.removeFromSubset(c, u);
		}
	}
}

Cover EgoSplitting::getCover() {
	return resultCover;
}

std::string EgoSplitting::toString() const {
	return "egosplitting";
}

std::unordered_map<node, std::vector<WeightedEdge>> EgoSplitting::getEgoNets() {
	return egoNets;
}

std::vector<std::unordered_map<node, index>> EgoSplitting::getEgoNetPartitions() {
	return egoNetExtendedPartitions;
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

void EgoSplitting::cleanUpCover() {
	if (parameters.at("Cleanup") == "Yes") {
		discardSmallCommunities();
		SignificanceCommunityCleanUp cleanup(G, resultCover);
		cleanup.run();
		resultCover = cleanup.getCover();
	}
	discardSmallCommunities();
}

} /* namespace NetworKit */
