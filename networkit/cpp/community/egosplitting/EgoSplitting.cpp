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

#include <tlx/unused.hpp>

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
#include "../../auxiliary/Parallelism.h"

namespace NetworKit {

EgoSplitting::EgoSplitting(const Graph &G, bool parallelEgoNetEvaluation)
		: EgoSplitting(G, parallelEgoNetEvaluation,
		               PLMFactory(true, 1.0, "none randomized").getFunction(),
		               LouvainMapEquationFactory(true, 16, "RelaxMap").getFunction()
) {
	INFO("Default EgoSplitting");
}

EgoSplitting::EgoSplitting(const Graph &G, bool parallelEgoNetEvaluation, ClusteringFunction clusterAlgo)
		: EgoSplitting(G, parallelEgoNetEvaluation, clusterAlgo, std::move(clusterAlgo)) {
}

EgoSplitting::EgoSplitting(const Graph &G, bool parallelEgoNetEvaluation, ClusteringFunction localClusterAlgo,
                           ClusteringFunction globalClusterAlgo)
		: Algorithm(), G(G),
		  parallelEgoNetEvaluation(parallelEgoNetEvaluation),
		  localClusteringAlgo(std::move(localClusterAlgo)),
		  globalClusteringAlgo(std::move(globalClusterAlgo)),
		  stochasticDistribution(0) {
	init();
}

void EgoSplitting::init() {
	personaEdges.resize(G.upperNodeIdBound());
	egoNetPartitions.resize(G.upperNodeIdBound());
	egoNetPartitionCounts.resize(G.upperNodeIdBound(), 0);
	personaOffsets.resize(G.upperNodeIdBound() + 1, 0);
	directedG = LowToHighDirectedGraph(G);

	parameters["storeEgoNet"] = "No";
	parameters["partitionFromGroundTruth"] = "No";
	parameters["numEgoNetsStored"] = "2000";
	parameters["Cleanup"] = "Yes";
	parameters["CleanupMerge"] = "Yes";
	parameters["maxEgoNetsPartitioned"] = "-1";

	// Connect Personas
	parameters["connectPersonas"] = "Yes";
	parameters["normalizePersonaCut"] = "No";
	parameters["connectPersonasStrat"] = "spanning";
	parameters["normalizePersonaWeights"] = "unweighted";

	// Parameters for ego-net extension
	parameters["Extend EgoNet Strategy"] = "Edges";
	parameters["Maximum Extend Factor"] = "5";
	parameters["addNodesExponent"] = "0.5";
	parameters["minNodeDegree"] = "2";
	parameters["Extend and Partition Iterations"] = "1";

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

	if (parameters.at("Cleanup") == "Yes"
	    || (parameters.at("Extend EgoNet Strategy") == "Edges" && parameters.at("Edges Score Strategy") == "Significance")
	    || parameters.at("Extend EgoNet Strategy") == "Significance") {
		stochasticDistribution.increaseMaxValueTo(2 * G.numberOfEdges() + G.numberOfNodes());
	}
	
	if (parameters.at("storeEgoNet") == "Yes") {
		egoNetExtendedPartitions.resize(G.upperNodeIdBound());
	}

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

	INFO("create Communities");
	std::vector<std::vector<node>> communities = getCommunitiesFromPersonaClustering();
	personaGraph = Graph(); // deallocate persona graph
	addTime(timer, "5  Create Communities");

	INFO("clean up Communities");
	cleanUpCommunities(communities);
	addTime(timer, "6  Cleanup Communities");

	INFO("create Cover");
	createCover(communities);
	addTime(timer, "7  Create Cover");

	hasRun = true;
}


void EgoSplitting::createEgoNets() {
	const int maxEgoNetsPartitioned = std::stoi(parameters.at("maxEgoNetsPartitioned"));
#pragma omp parallel if (parallelEgoNetEvaluation)
	{
		Aux::SignalHandler signalHandler;
		//Aux::Timer timer, totalTimer;
		//timer.start();
		//totalTimer.start();
		NodeMapping egoMapping(G); // Assign local IDs to the neighbors
		MemoizationTable<double> sigTable(-1.0, 0);
		SparseVector<double> nodeScores(0);
		SparseVector<node> significantGroup(0, none);
		SparseVector<std::vector<count>> edgesToGroups(0);
		SignificanceCalculator significance(stochasticDistribution);
		EgoNetData egoNetData{G, directedG, groundTruth, parameters, sigTable, egoMapping, nodeScores,
			significantGroup, edgesToGroups, significance};
		//addTime(timer, "11    Data Setup");

#pragma omp for schedule(dynamic, 10)
		for (omp_index egoNode = 0; egoNode < static_cast<omp_index>(G.upperNodeIdBound()); ++egoNode) {
			if (!signalHandler.isRunning()) continue;

			if (!G.hasNode(egoNode) || G.degree(egoNode) < 2)
				continue;

			if (maxEgoNetsPartitioned != -1 && egoNode > maxEgoNetsPartitioned) {
				egoNetPartitionCounts[egoNode] = 1;
				G.forNeighborsOf(egoNode, [&](node neighbor) {
					egoNetPartitions[egoNode].emplace(neighbor, 0);
				});
				continue;
			}

			egoMapping.reset();

			DEBUG("Create EgoNet for Node ", egoNode, "/", G.upperNodeIdBound());
			// Find neighbors == nodes of the ego-net
			G.forEdgesOf(egoNode, [&](node, node v) {
				if (G.degree(v) > 1) {
					egoMapping.addNode(v);
				}
			});

			const count originalEgoNetSize = egoMapping.nodeCount();
			tlx::unused(originalEgoNetSize);

			Graph egoGraph(egoMapping.nodeCount(), G.isWeighted());

			// Find all triangles and add the edges to the egoGraph
			egoGraph.forNodes([&](const node localV) {
				const node globalV = egoMapping.toGlobal(localV);
				directedG.forEdgesOf(globalV, [&](node, node w, edgeweight weight2) {
					if (egoMapping.isMapped(w)) {
						// we have found a triangle u-v-w
						egoGraph.addEdge(localV, egoMapping.toLocal(w), weight2);
					}
				});
			});

			// Handle some more trivial cases: egonet is a clique or just two unconnected nodes
			{
				// egoNode has only degree-1-neighbors
				if (egoGraph.numberOfNodes() == 0) continue;

				// Test if the egonet is a clique, if yes, do not split the egonet
				count cliqueEdges = egoGraph.numberOfNodes() * (egoGraph.numberOfNodes() - 1) / 2;
				if (egoGraph.numberOfEdges() == cliqueEdges) {
					egoNetPartitionCounts[egoNode] = 1;
					for (node neighbor : egoMapping.globalNodes()) {
						egoNetPartitions[egoNode].emplace(neighbor, 0);
					}
					continue;
				}
			}

			//addTime(timer, "13    Build EgoNet");

			DEBUG("Extend and partition");
			EgoNetExtensionAndPartition extAndPartition(egoNetData, egoNode, egoGraph,
			                                            localClusteringAlgo);
			extAndPartition.run();
			const Partition &egoPartition = extAndPartition.getPartition();
			//addTimings(extAndPartition.getTimings(), "15");
			//addTime(timer, "15    Extend and Partition EgoNet");

			if (parameters.at("storeEgoNet") == "Yes") { // only for analysis
				Graph extendedEgoGraph = extAndPartition.getExtendedEgoGraph();
				storeEgoNet(extendedEgoGraph, egoMapping, egoNode);
				// Store ego-net partition with extended nodes
				extendedEgoGraph.forNodes([&](node i) {
					egoNetExtendedPartitions[egoNode].emplace(egoMapping.toGlobal(i),
					                                          egoPartition.subsetOf(i));
				});
			}
			//addTime(timer, "18    Store EgoNet");

			DEBUG("Store ego-net");
			// Remove nodes that are not directed neighbors (they were added by the ego-net extension)
			Partition directNeighborPartition(egoGraph.upperNodeIdBound());
			directNeighborPartition.setUpperBound(egoPartition.upperBound());
			egoGraph.forNodes([&](node v) {
				directNeighborPartition.addToSubset(egoPartition.subsetOf(v), v);
			});
			directNeighborPartition.compact(true);
			egoNetPartitionCounts[egoNode] = directNeighborPartition.upperBound();
			assert(egoNetPartitionCounts[egoNode] == directNeighborPartition.numberOfSubsets());
			egoGraph.forNodes([&](const node localI) {
				const node i = egoMapping.toGlobal(localI);
				egoNetPartitions[egoNode].emplace(i, directNeighborPartition[localI]);
			});

			assert(G.degree(egoNode) >= originalEgoNetSize);
			assert(egoNetPartitions[egoNode].size() == originalEgoNetSize);
			assert(egoGraph.numberOfNodes() == originalEgoNetSize);
			//addTime(timer, "19    Store EgoNet Partition");

			if (parameters.at("connectPersonas") == "Yes")
				personaEdges[egoNode] = connectEgoPartitionPersonas(egoGraph, directNeighborPartition);
			//addTime(timer, "16    Connect Ego-Personas");

			//addTime(timer, "1x    Clean up");
		}
		//addTime(totalTimer, "10    EgoNet Sum");
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
	DEBUG("Connect personas of one node");
	std::vector<WeightedEdge> edges;

	// Contract graph
	ParallelPartitionCoarsening coarsening{egoGraph, egoPartition, false, false};
	coarsening.run();
	const Graph &coarseGraph = coarsening.getCoarseGraph();
	const std::vector<node> &egoToCoarse = coarsening.getFineToCoarseNodeMapping();

	// Build a mapping from nodes in the coarse graph to partition ids in egoPartition
	// Note that if we could assume that calling Partition::compact() twice does not modify partition ids, we might be able to omit this altogether
	std::vector<node> personaIndex(coarseGraph.upperNodeIdBound(), none);
	egoGraph.forNodes([&](node u) {
		personaIndex[egoToCoarse[u]] = egoPartition[u];
	});

	auto addPersonaEdge = [&](node u, node v, edgeweight w = 1) {
		index p_u = personaIndex[u];
		index p_v = personaIndex[v];
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
		const Graph &spanningForest = span.getMSF();
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
		sum += egoNetPartitionCounts[u];
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
	// Connect personas of different nodes
	count twiceNumberOfEdges = 0;
	#pragma omp parallel
	{
		count localTwiceNumberOfEdges = 0;

		#pragma omp for schedule (dynamic, 1000) nowait
		for (node u = 0; u < G.upperNodeIdBound(); ++u) {
			if (G.hasNode(u)) {
				if (G.degree(u) >= 2) {
					for (const WeightedEdge& edge : personaEdges[u]) {
						node pu = getPersona(u, edge.u);
						node pv = getPersona(u, edge.v);
						assert(pu != pv);
						personaGraph.addHalfEdge(pu, pv, edge.weight);
						personaGraph.addHalfEdge(pv, pu, edge.weight);
						localTwiceNumberOfEdges += 2;
					}

					G.forEdgesOf(u, [&](node , node v, edgeweight weight) {
						if (G.degree(v) >= 2) {
							auto idx_u = egoNetPartitions[u].find(v);
							auto idx_v = egoNetPartitions[v].find(u);
							assert(idx_u != egoNetPartitions[u].end() && idx_v != egoNetPartitions[v].end());
							const node pu = getPersona(u, idx_u->second);
							const node pv = getPersona(v, idx_v->second);
							assert(pu != pv);
							personaGraph.addHalfEdge(pu, pv, weight);
							++localTwiceNumberOfEdges;
						}
					});
				}
			}
		}
		
		#pragma omp atomic update
		twiceNumberOfEdges += localTwiceNumberOfEdges;
	}
	
	if (twiceNumberOfEdges % 2 != 0)
		throw std::runtime_error("edge count accumulation failed");
	personaGraph.setNumberOfEdges(twiceNumberOfEdges / 2);
	
#ifndef NDEBUG
	count internalPersonaEdges = 0;
	for (const auto &edges : personaEdges)
		internalPersonaEdges += edges.size();
	count removedEdges = 0;
	G.forNodes([&](node v) {
		if (G.degree(v) == 1) {
			node neighbor = none;
			G.forNeighborsOf(v, [&](node x) { neighbor = x; });
			assert(neighbor != none);
			if (G.degree(neighbor) > 1 || v < neighbor) {
				++removedEdges;
			}
		}
	});
	assert(personaGraph.numberOfEdges() + removedEdges == G.numberOfEdges() + internalPersonaEdges);

	personaGraph.forNodes([&](node u) {
		assert(personaGraph.degree(u) >= 1);
		if (personaGraph.degree(u) == 1) {
			node neighbor = none;
			personaGraph.forNeighborsOf(u, [&](node x) { neighbor = x; });
			assert(neighbor != u);
		}
	});
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
	auto iso2 = numIsolatedNodes(personaGraph);
	assert(iso2 == 0);
#endif

	ConnectedComponents compsAlgo(personaGraph);
	compsAlgo.run();
	std::vector<count> componentSizes(compsAlgo.numberOfComponents());
	personaGraph.forNodes([&](node u) {
		++componentSizes[compsAlgo.componentOfNode(u)];
	});

	personaGraph.forNodes([&](node u) {
		const index c = compsAlgo.componentOfNode(u);
		if (componentSizes[c] < 5) {
			personaGraph.removeNode(u);
		}
	});
	
	egoNetPartitions.clear();
	egoNetPartitions.shrink_to_fit();
	
}

void EgoSplitting::createPersonaClustering() {
	personaPartition = globalClusteringAlgo(personaGraph);
}


std::vector<std::vector<node>> EgoSplitting::getCommunitiesFromPersonaClustering() {
	personaPartition.compact(personaPartition.upperBound() < 2 * personaGraph.numberOfNodes());
	std::vector<std::vector<node>> result(personaPartition.upperBound());
	G.forNodes([&](node u) {
		for (index i = personaOffsets[u]; i < personaOffsets[u + 1]; ++i) {
			if (personaGraph.hasNode(i)) {
				index part = personaPartition.subsetOf(i);
				if (result[part].empty() || result[part].back() != u) {
					result[part].push_back(u);
				}
			}
		}
	});

	return result;
}

void EgoSplitting::createCover(const std::vector<std::vector<node>> &communities) {
	// Create cover from given communities
	resultCover = Cover(G.upperNodeIdBound());
	resultCover.setUpperBound(communities.size());
	for (index com_id = 0; com_id < communities.size(); ++com_id) {
		for (node u : communities[com_id]) {
			resultCover.addToSubset(com_id, u);
		}
	}
}

void EgoSplitting::cleanUpCommunities(std::vector<std::vector<node>> &communities) {
	if (parameters.at("Cleanup") == "Yes") {
		discardSmallCommunities(communities);
		bool mergeDiscarded = parameters.at("CleanupMerge") == "Yes";
		SignificanceCommunityCleanUp cleanup(G, communities, stochasticDistribution, 0.1, 0.1, 0.5, mergeDiscarded);
		cleanup.run();
	}
	discardSmallCommunities(communities);
}

void EgoSplitting::discardSmallCommunities(std::vector<std::vector<node>> &communities) {// Discard communities of size 4 or less
	count min_size = 5;

	auto new_end = std::remove_if(communities.begin(), communities.end(),
				      [min_size](const std::vector<node> &c) {
					      return c.size() < min_size;
				      });
	communities.erase(new_end, communities.end());
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

} /* namespace NetworKit */
