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

#include "../structures/UnionFind.h"
#include "../structures/Partition.h"
#include "../components/ConnectedComponents.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/SignalHandling.h"
#include "../auxiliary/Timer.h"
#include "../sparsification/LocalSimilarityScore.h"
#include "../edgescores/TriangleEdgeScore.h"
#include "../edgescores/EdgeScoreAsWeight.h"
#include "../coarsening/ParallelPartitionCoarsening.h"
#include "../graph/RandomMaximumSpanningForest.h"
#include "EgoSplitting.h"
#include "PLM.h"
#include "CoverF1Similarity.h"

#define true_or_throw(cond, msg) if (!cond) throw std::runtime_error(msg)

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
	directedEdges = AdjacencyArray(G);

	parameters["storeEgoNet"] = "No";
	parameters["addEgoNode"] = "No";
	parameters["weightedEgoNet"] = "No";
	parameters["partitionFromGroundTruth"] = "No";
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
	parameters["processEgoNet"] = "extend";
	parameters["addNodesFactor"] = "4";
	parameters["addNodesExponent"] = "0.6";
	parameters["edgesBetweenNeigNeig"] = "Yes";
	parameters["keepOnlyTriangles"] = "No";
	parameters["minNodeDegree"] = "0";
	parameters["triangleThreshold"] = "0";

	// Parameters for edgeScores
	parameters["scoreStrategy"] = "count";
	parameters["extendRandom"] = "No";
	parameters["extendOverDirected"] = "No";

	// Test parameters
	parameters["processEgoNet"] = "none";
//	parameters["storeEgoNet"] = "Yes";
//	parameters["processEgoNet"] = "extend";
//	parameters["extendRandom"] = "No";
//	parameters["extendStrategy"] = "edgeScore";
//	parameters["scoreStrategy"] = "score^2_normed";
//	parameters["minTriangles"] = "2";
//	parameters["storeEgoNet"] = "Yes";
//	parameters["triangleThreshold"] = "0";
//	parameters["minNodeDegree"] = "2";
//	parameters["partitionFromGroundTruth"] = "Yes";
}

void EgoSplitting::run() {
	Aux::SignalHandler handler;
	Aux::Timer timer;
	timings.clear();

	if (std::stod("0.3") == 0)
		throw std::runtime_error("Can't convert numbers because of wrong locale!");
	G.forNodes([&](node u) {
		if (G.hasEdge(u, u))
			throw std::runtime_error("No self-loops allowed!");
	});


//	edgeScoreGraph = weightedEdgesGraph(G);

	INFO("create EgoNets");
	timer.start();
	createEgoNets();
	timer.stop();
	timings["1) create EgoNets"] = timer.elapsedMicroseconds();
	handler.assureRunning();

	INFO("split into Personas");
	timer.start();
	splitIntoPersonas();
	timer.stop();
	timings["2) split Personas"] = timer.elapsedMicroseconds();
	handler.assureRunning();

	INFO("connect Personas");
	timer.start();
	connectPersonas();
	timer.stop();
	timings["3) connect Personas"] = timer.elapsedMicroseconds();
	handler.assureRunning();

	INFO("create Persona Clustering");
	timer.start();
	createPersonaClustering();
	timer.stop();
	timings["4) Persona Clustering"] = timer.elapsedMicroseconds();
	handler.assureRunning();

	INFO("create Cover");
	timer.start();
	createCover();
	timer.stop();
	timings["5) create Cover"] = timer.elapsedMicroseconds();

	hasRun = true;
}

std::string EgoSplitting::toString() const {
	return "EgoSplitting";
}

Graph EgoSplitting::weightedEdgesGraph(Graph const &inputGraph) {
	// Get triangle scores
	auto trianglesAlgo = TriangleEdgeScore(inputGraph);
	trianglesAlgo.run();
	auto triangleScores = trianglesAlgo.scores();
	// Calculate similarity scores
	auto localSimScoreAlgo = LocalSimilarityScore(inputGraph, triangleScores);
	localSimScoreAlgo.run();
	auto localSimScores = localSimScoreAlgo.scores();
	// Create weighted graph from scores
	double weightFactor = std::stod(parameters.at("weightFactor"));
	double weightOffset = std::stod(parameters.at("weightOffset"));
	auto toWeight = EdgeScoreAsWeight(inputGraph, localSimScores, false, weightOffset,
	                                  weightFactor);
	Graph weightedGraph = toWeight.calculate();
	return weightedGraph;
}

void EgoSplitting::createEgoNets() {
	// Assign IDs to the neighbours
	NodeMapping nodeMapping(G);

	Aux::SignalHandler handler;
	Aux::Timer timer;

	G.forNodes([&](node u) {
		handler.assureRunning();
		INFO("Create EgoNet for Node ", u, "/", G.upperNodeIdBound());

		/******************************************************************************************
		 **                                Add Neighbor Nodes                                    **
		 ******************************************************************************************/
		timer.start();
		count degree = G.degree(u);
		Graph egoGraph(degree, true);

		if (parameters.at("addEgoNode") == "Yes") {
			nodeMapping.addNode(u);
			egoGraph.addNode();
		}
		// Add neighbors
		G.forEdgesOf(u, [&](node, node v) {
			nodeMapping.addNode(v);
		});
		timer.stop();
		timings["1a)    Find nodes"] += timer.elapsedMicroseconds();


		/******************************************************************************************
		 **                             Triangle Search for Edges                                **
		 ******************************************************************************************/
		timer.start();
		// Find all triangles and add the edges to the egoGraph
		G.forEdgesOf(u, [&](node, node v, edgeweight weight1) {
			if (parameters.at("addEgoNode") == "Yes") {
				egoGraph.addEdge(nodeMapping.local(u), nodeMapping.local(v), weight1);
			}
			directedEdges.forEdgesOf(v, [&](node, node w, edgeweight weight2) {
				if (nodeMapping.isMapped(w)) {
					// we have found a triangle u-v-w
					egoGraph.addEdge(nodeMapping.local(v), nodeMapping.local(w), weight2);
				}
			});
		});
		timer.stop();
		timings["1b)    Neighbor Triangle Search"] += timer.elapsedMicroseconds();


		/******************************************************************************************
		 **                              Extend/Process EgoNet                                   **
		 ******************************************************************************************/
		timer.start();
		double addNodesFactor = std::stod(parameters.at("addNodesFactor"));
		double addNodesExponent = std::stod(parameters.at("addNodesExponent"));
		count extendNodeCnt = std::ceil(
				addNodesFactor * std::pow(degree, addNodesExponent));
		if (parameters.at("processEgoNet") == "none") {
			// nothing to do
		} else if (parameters.at("processEgoNet") == "extend") {
			extendEgoNet(egoGraph, u, nodeMapping, extendNodeCnt);
		} else if (parameters.at("processEgoNet") == "consensus") {
			consensusWeighting(egoGraph, u, nodeMapping, extendNodeCnt);
		} else {
			throw std::runtime_error("Missing strategy to extend ego-net");
		}
		// TODO: score basierend auf edges / degree des Kandidaten
		// TODO: Überlappende Partitionierung
		timer.stop();
		timings["1c)    Extend EgoNet"] += timer.elapsedMicroseconds();

		/******************************************************************************************
		 **                                  Cluster EgoNet                                      **
		 ******************************************************************************************/
		timer.start();
		// Cluster ego-net with the local cluster algorithm
		Partition egoPartition;
		if (parameters.at("weightedEgoNet") != "Yes")
			egoGraph = Graph(egoGraph, false, egoGraph.isDirected());
		if (parameters.at("partitionFromGroundTruth") == "Yes")
			egoPartition = createGroundTruthPartition(egoGraph, nodeMapping, u);
		else if (egoGraph.numberOfEdges() > 0) {
			bool weighted = parameters.at("weightedEgoNet") == "Yes";
//			Graph egoCopy(egoGraph, weighted, egoGraph.isDirected());
			egoPartition = localClusterAlgo(egoGraph); // Python moves the graph so we need a copy
		} else {
			egoPartition = Partition(degree);
			egoPartition.allToSingletons();
		}
		egoPartition.compact();
		timer.stop();
		timings["1d)    Cluster EgoNet"] += timer.elapsedMicroseconds();


		/******************************************************************************************
		 **                                 Connect Personas                                     **
		 ******************************************************************************************/
		if (parameters.at("connectPersonas") == "Yes")
			personaEdges[u] = connectEgoPartitionPersonas(egoGraph, egoPartition);



		/******************************************************************************************
		 **                                  Store EgoNet                                        **
		 ******************************************************************************************/
		timer.start();
		if (parameters.at("storeEgoNet") == "Yes") { // && u < 100
			// Get EgoNet with gloabl node ids
			Graph egoNetGraph = Graph(G.upperNodeIdBound(), egoGraph.isWeighted());
			egoGraph.forEdges([&](node v, node w, edgeweight weight) {
				egoNetGraph.addEdge(nodeMapping.global(v), nodeMapping.global(w), weight);
			});
			for (node v : egoNetGraph.nodes()) {
				if (!egoGraph.hasNode(nodeMapping.local(v)))
					egoNetGraph.removeNode(v);
			}
			if (egoGraph.isWeighted())
				egoNetGraph.addEdge(egoNets[u].nodes()[0], egoNets[u].nodes()[1],
				                   0.000001);
			egoNets[u] = egoNetGraph;
		}
		timer.stop();
		timings["1x)    Copy EgoNet"] += timer.elapsedMicroseconds();


		/******************************************************************************************
		 **                                 Build EgoNet Map                                     **
		 ******************************************************************************************/
		timer.start();
		// Insert nodes into ego-net data structure
		for (node i : egoGraph.nodes()) {
			egoNetPartitions[u].emplace(nodeMapping.global(i), egoPartition.subsetOf(i));
		}
		assert(egoNetPartitions[u].size() == egoGraph.numberOfNodes());
		egoNetPartitionCounts[u] = egoPartition.numberOfSubsets();
		timer.stop();
		timings["1e)    EgoNet subsets"] += timer.elapsedMicroseconds();


		/******************************************************************************************
		 **                                     Cleanup                                          **
		 ******************************************************************************************/
		timer.start();
		nodeMapping.reset();
		timer.stop();
		timings["1f)    Clean up"] += timer.elapsedMicroseconds();
	});
}

std::vector<std::tuple<index, index, edgeweight>>
EgoSplitting::connectEgoPartitionPersonas(const Graph &egoGraph, const Partition &egoPartition) {
	std::vector<std::tuple<index, index, edgeweight >> edges;
	// Contract graph
	ParallelPartitionCoarsening coarsening(egoGraph, egoPartition);
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
	double spanSize;
	{
		// Get spanning forest
		RandomMaximumSpanningForest span(coarseGraph);
		span.run();
		auto spanningForest = span.getMSF();
		spanSize = spanningForest.numberOfEdges();
	}

	if (parameters.at("normalizePersonaCut") == "volume") {
		coarseGraph.forEdges([&](node u, node v, edgeweight w) {
			double volume = coarseGraph.weight(u, u) + coarseGraph.weight(v, v);
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

	std::basic_string<char> &strategy = parameters.at("connectPersonasStrat");
	if (strategy == "spanning") {
		// TODO: Better weights
		// Every persona gets at most 'iterations' edges
		count iterations = std::stoi(parameters.at("maxPersonaEdges"));
		for (int i = 0; i < iterations; ++i) {
			RandomMaximumSpanningForest span(coarseGraph);
			span.run();
			auto spanningForest = span.getMSF();
			spanningForest.forEdges([&](node u, node v, edgeweight w) {
				if (parameters.at("iterationWeight") == "Yes")
					w = w / (i + 1);
				addPersonaEdge(u, v, w);
				coarseGraph.removeEdge(u, v);
			});
		}
	} else if (strategy == "maxEdge") {
		count iterations = std::stoi(parameters.at("maxPersonaEdges"));
		for (int i = 0; i < iterations; ++i) {
			std::vector<std::set<node>> edgeSets(coarseGraph.upperNodeIdBound());
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
		throw std::runtime_error(strategy + " is no valid strategy to connect personas!");
	}

	if (parameters.at("normalizePersonaWeights") == "spanSize") {
		edgeweight weightSum = 0.0;
		for (auto &edge : edges)
			weightSum += std::get<2>(edge);
		weightSum /= spanSize; // Sum of edge weights == spanning tree size
		for (auto &edge : edges) {
			std::get<2>(edge) /= weightSum;
		}
	} else if (parameters.at("normalizePersonaWeights") == "unweighted") {
		for (auto &edge : edges)
			std::get<2>(edge) = 1;
	} else if (parameters.at("normalizePersonaWeights") == "sameWeights") {
		for (auto &edge : edges)
			std::get<2>(edge) = spanSize / edges.size();
	}
	return edges;
}

void
EgoSplitting::consensusWeighting(Graph &egoGraph, node u, const NodeMapping &nodeMapping,
                                 count extendNodeCnt) {
	if (parameters.at("extendRandom") != "Yes")
		throw std::runtime_error("Consensus is only useful with randomization");
	Graph egoGraphWeighted(egoGraph, true, egoGraph.isDirected());
	egoGraphWeighted.forEdges([&](node v, node w) {
		egoGraphWeighted.setWeight(v, w, 0);
	});
	int iterations = 10;
	for (int it = 0; it < iterations; ++it) {
		Graph extendedEgoGraph = egoGraph;
		NodeMapping extendedNodeMapping = nodeMapping;
		extendEgoNet(extendedEgoGraph, u, extendedNodeMapping,
		             extendNodeCnt);
		auto egoPartition = localClusterAlgo(extendedEgoGraph);
		egoGraphWeighted.forEdges([&](node v, node w, edgeweight weight) {
			if (egoPartition.inSameSubset(v, w))
				egoGraphWeighted.increaseWeight(v, w, 1);
		});
	}
	egoGraph = egoGraphWeighted;
}

void EgoSplitting::findTriangles(Graph graph, AdjacencyArray directedGraph,
                                 std::function<void(node, node,
                                                    node)> triangleFunc) const {
	std::vector<int> currentNeighbors(graph.upperNodeIdBound(), 0);
	graph.forNodes([&](node v) {
		std::vector<node> currentNeighborsList;
		graph.forEdgesOf(v, [&](node, node w, edgeweight weight1) {
			currentNeighbors[w] = 1;
			currentNeighborsList.push_back(w);
		});
		directedGraph.forEdgesOf(v, [&](node, node w, edgeweight weight1) {
			directedGraph.forEdgesOf(w, [&](node, node x, edgeweight weight2) {
				if (currentNeighbors[x] == 0)
					return;
				// We found a triangle v-w-x
				triangleFunc(v, w, x);
			});
		});
		for (node w : currentNeighborsList)
			currentNeighbors[w] = 0;
	});
}

void
EgoSplitting::extendEgoNet(Graph &egoGraph, node u, NodeMapping &neighbors,
                           count extendNodeCnt) {
	const count directNeighborsCnt = neighbors.nodeCount();
	std::vector<std::pair<node, double>> nodeScores; // node and its score

	/**********************************************************************************************
	 **                           Get node candidates with scores                                **
	 **********************************************************************************************/
	const auto extendStrategy = parameters.at("extendStrategy");

	std::vector<std::set<node>> triangleEdges;
	if (extendStrategy == "edgeScore")
		nodeScores = scoreEdgeCount(u, neighbors);
	else if (extendStrategy == "triangles")
		nodeScores = scoreTriangles(u, neighbors, triangleEdges, egoGraph);
	else
		throw std::runtime_error("No valid strategy to extend Ego-Net");

	// Remove nodes with score zero
	for (size_t i = 0; i < nodeScores.size(); ++i) {
		if (nodeScores[i].second <= 0.0) {
			nodeScores[i] = nodeScores.back();
			nodeScores.pop_back();
			--i; // Check node at current position again
		}
	}
	// TODO: Use only nodes with score > 0.1 (or something), use all nodes with score > 0.5
	// TODO: Score needs to be normalized to (0, 1) for this to work


	/**********************************************************************************************
	 **                                  Add nodes to result                                     **
	 **********************************************************************************************/
	std::vector<node> nodesToAdd;
	if (parameters.at("extendRandom") == "Yes") {
		// Take nodes randomly, the chance to pick a node is proportional to its score
		std::vector<double> random_weights;
		for (auto pair : nodeScores)
			random_weights.push_back(pair.second);

		for (index i = 0; i < extendNodeCnt && !nodeScores.empty(); ++i) {
			std::discrete_distribution<index> dist(random_weights.begin(),
			                                       random_weights.end());
			index random_node = dist(Aux::Random::getURNG());
			nodesToAdd.push_back(random_node);
			nodeScores[random_node] = nodeScores.back();
			nodeScores.pop_back();
			random_weights[random_node] = random_weights.back();
			random_weights.pop_back();
		}
	} else {
		// Take the nodes with the best scores
		std::sort(nodeScores.begin(), nodeScores.end(),
		          [](std::pair<node, double> a, std::pair<node, double> b) {
			          return a.second > b.second;
		          });
		if (nodeScores.size() > extendNodeCnt)
			nodeScores.resize(extendNodeCnt);
		for (auto pair : nodeScores) {
			nodesToAdd.push_back(pair.first);
		}
	}
	for (node v : nodesToAdd) {
		neighbors.addNode(v);
		egoGraph.addNode();
	}
	Graph onlyNeighborsGraph(egoGraph);


	/**********************************************************************************************
	 **                                  Add edges to result                                     **
	 **********************************************************************************************/
	for (node v : neighbors.globalNodes()) {
		directedEdges.forEdgesOf(v, [&](node, node w, edgeweight weight) {
			if (neighbors.isMapped(w)) {
				node v_loc = neighbors.local(v);
				node w_loc = neighbors.local(w);
				// Edges between direct neighbors are already in the Egonet
				if (v_loc < directNeighborsCnt && w_loc < directNeighborsCnt)
					return;

				// Discard edges between neighbors of neighbors
				if (parameters.at("edgesBetweenNeigNeig") != "Yes"
				    && v_loc >= directNeighborsCnt
				    && w_loc >= directNeighborsCnt) {
					return;
				}
				egoGraph.addEdge(v_loc, w_loc, weight);
			}
		});
	}


	/**********************************************************************************************
	 **                         Use only triangle edges (optional)                               **
	 **********************************************************************************************/
	if (parameters.at("keepOnlyTriangles") == "Yes") {
		AdjacencyArray directedEgoGraph(egoGraph);
		auto foundTriangleWork = [&](node v, node w, node x) {
			node v_loc = neighbors.local(v);
			node w_loc = neighbors.local(w);
			node x_loc = neighbors.local(x);
			if (v_loc >= directNeighborsCnt || w_loc >= directNeighborsCnt)
				onlyNeighborsGraph.addEdge(v, w);
			if (w_loc >= directNeighborsCnt || x_loc >= directNeighborsCnt)
				onlyNeighborsGraph.addEdge(w, x);
			if (x_loc >= directNeighborsCnt || v_loc >= directNeighborsCnt)
				onlyNeighborsGraph.addEdge(x, v);
		};
		findTriangles(egoGraph, directedEgoGraph, foundTriangleWork);
		egoGraph = onlyNeighborsGraph;
//		for (node v : nodesToAdd) {
//			for (node w : triangleEdges[v])
//				egoGraph.addEdge(neighbors.local(v), neighbors.local(w), 1);
//		}
	}


	/**********************************************************************************************
	 **           Weight graph based on number of triangles for each edge (optional)             **
	 **********************************************************************************************/
	 if (parameters.at("weightedEgoNet") == "Yes")	{
		Graph egoWeightsGraph(egoGraph, true, false);
		egoWeightsGraph.forEdges([&](node v, node w, edgeweight weight) {
			egoWeightsGraph.setWeight(v, w, 0);
		});
		AdjacencyArray directedEgoGraph(egoGraph);
		auto foundTriangleWork = [&](node v, node w, node x) {
			egoWeightsGraph.increaseWeight(v, w, 1);
			egoWeightsGraph.increaseWeight(w, x, 1);
			egoWeightsGraph.increaseWeight(x, v, 1);
		};
		findTriangles(egoGraph, directedEgoGraph, foundTriangleWork);

		egoGraph = egoWeightsGraph;
	}

	count minDegree = stoi(parameters.at("minNodeDegree"));
	removeLowDegreeNodes(egoGraph, minDegree, directNeighborsCnt);
	removeLowTriangleCntNodes(egoGraph, directNeighborsCnt); // TODO: Bevorzuge high degree nodes (?)
	removeLowDegreeNodes(egoGraph, minDegree, directNeighborsCnt);
}


void
EgoSplitting::removeLowTriangleCntNodes(Graph &egoGraph, count directNeighborsCnt) const {
	double triangleThreshold = std::stod(parameters.at("triangleThreshold"));
	if (triangleThreshold <= 0)
		return;
	AdjacencyArray directedEgoGraph(egoGraph);
	count nodes_changed = 1;
	std::vector<count> triangleCount(egoGraph.upperNodeIdBound(), 0);
	std::vector<int> currentNeighbors(G.upperNodeIdBound(), 0);
	auto foundTriangleWork = [&](node v, node w, node x) {
		++triangleCount[v];
		++triangleCount[w];
		++triangleCount[x];
	};
	while (nodes_changed > 0) {
		nodes_changed = 0;
		findTriangles(egoGraph, directedEgoGraph, foundTriangleWork);
		egoGraph.forNodesInRandomOrder([&](node v) {
			if (v >= directNeighborsCnt && triangleCount[v] * 1.0 / egoGraph.degree(v) < triangleThreshold) {
				++nodes_changed;
				egoGraph.removeNode(v);
			}
		});
	}
}

void EgoSplitting::removeLowDegreeNodes(Graph &egoGraph, count minDegree,
                                        count directNeighborsCnt) const {
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

double EgoSplitting::normalizeScore(node v, double score) const {
	std::string scoreStrategy = parameters.at("scoreStrategy");
	if (scoreStrategy == "constant")
		return 1.0;
	if (scoreStrategy == "score")
		return score * 1.0;
	if (scoreStrategy == "score_normed")
		return score * 1.0 / G.degree(v);
	if (scoreStrategy == "score^1.5_normed")
		return std::pow(score, 1.5) / G.degree(v);
	if (scoreStrategy == "score^2_normed")
		return score * score * 1.0 / G.degree(v);
	throw std::runtime_error("No valid score strategy provided");
}

std::vector<std::pair<node, double>>
EgoSplitting::scoreEdgeCount(node u, const NodeMapping &neighbors) {
	std::vector<std::vector<double>> edgeScores(G.upperNodeIdBound());
	std::vector<node> secondNeighbors;
	// Search for all edges to neighbors of neighbors
	std::vector<node> directNeighbors = neighbors.globalNodes();
	auto isDirectNeighbor = [&](node x) {
		return neighbors.isMapped(x);
	};
	auto work = [&](node v, node w, edgeweight weight) {
		if (!isDirectNeighbor(w) && w != u) {
			if (edgeScores[w].empty())
				secondNeighbors.push_back(w);
			edgeScores[w].push_back(weight);
		}
	};
	for (node v : directNeighbors) {
		if (parameters.at("extendOverDirected") == "Yes")
			directedEdges.forEdgesOf(v, work);
		else
			G.forEdgesOf(v, work);
	}
	// Calculate score for each candidate
	std::vector<std::pair<node, double>> nodeScores;
	for (node v : secondNeighbors) {
		double numEdges = edgeScores[v].size();
		double score = normalizeScore(v, numEdges);
		nodeScores.emplace_back(v, score);
	}
	return nodeScores;
}

std::vector<std::pair<node, double>>
EgoSplitting::scoreTriangles(node u, const NodeMapping &neighbors,
                             std::vector<std::set<node>> &triangleEdges,
                             Graph const &egoGraph) {
	std::vector<int> triangleCnts(G.upperNodeIdBound(), -1);
	std::vector<node> directNeighbors = neighbors.globalNodes();
	std::vector<node> secondNeighbors;
	auto isDirectNeighbor = [&](node x) {
		return neighbors.isMapped(x);
	};
	count directNeighborsCnt = directNeighbors.size();

	Graph extendedEgoGraph = egoGraph;
	NodeMapping extendedMapping = neighbors;

	for (node v : directNeighbors) {
		G.forEdgesOf(v, [&](node, node w, edgeweight weight) {
			if (!isDirectNeighbor(w) && w != u) {
				if (!extendedMapping.isMapped(w)) {
					secondNeighbors.push_back(w);
					triangleCnts[w] = 0;
					extendedEgoGraph.addNode();
					extendedMapping.addNode(w);
				}
				extendedEgoGraph.addEdge(extendedMapping.local(v),
				                         extendedMapping.local(w));
			}
		});
	}

	// If only one of the three nodes is not a direct neighbor, return that node
	auto getSingleSecondNeighborLocal = [&](node node1, node node2, node node3) {
		int b1 = (node1 >= directNeighborsCnt);
		int b2 = (node2 >= directNeighborsCnt);
		int b3 = (node3 >= directNeighborsCnt);
		if (b1 + b2 + b3 == 1) {
			if (b1)
				return node1;
			if (b2)
				return node2;
			if (b3)
				return node3;
		}
		return none;
	};

	AdjacencyArray directedEgoGraph(extendedEgoGraph);
	std::vector<int> currentNeighbors(G.upperNodeIdBound(), 0);
	triangleEdges.resize(G.upperNodeIdBound());
	auto foundTriangleWork = [&](node v, node w, node x) {
		// Only count triangles if exactly one node is a second neighbor
		node secondNeighbor = getSingleSecondNeighborLocal(v, w, x);
		if (secondNeighbor != none) {
			node globalNode = extendedMapping.global(secondNeighbor);
			++triangleCnts[globalNode];
			for (node neighborNode : {v, w, x}) {
				if (secondNeighbor != neighborNode) {
					triangleEdges[extendedMapping.global(secondNeighbor)].insert(
							extendedMapping.global(neighborNode));
				}
			}
		}
	};
	findTriangles(extendedEgoGraph, directedEgoGraph, foundTriangleWork);

	std::vector<std::pair<node, double>> nodeScores;
	nodeScores.reserve(secondNeighbors.size());
	for (node v : secondNeighbors) {
		double score = triangleCnts[v];
		// TODO: Maybe use number of triangle edges as score?
		// TODO: num_triangles / num_triangle_edges : Sinnvolle Metrik?
		// TODO: teile durch Anzahl Kanten ins Egonet (Anzahl mitzählen)
		// TODO: Gewichte Kanten mit Anzahl Triangles
		int minTriangles = std::stoi(parameters.at("minTriangles"));
		if (score < minTriangles)
			score = 0.0;
		double normalizedScore = normalizeScore(v, score);
		nodeScores.emplace_back(v, normalizedScore);
	}
	return nodeScores;
}

void EgoSplitting::splitIntoPersonas() {
	count sum = 0;
	for (index i = 0; i < G.upperNodeIdBound(); ++i) {
		personaOffsets[i] = sum;
		sum += egoNetPartitionCounts[i];
	}
	personaOffsets[G.upperNodeIdBound()] = sum;
	personaGraph = Graph(sum, true);

}

void EgoSplitting::connectPersonas() {
	auto getPersona = [&](node u, index i) {
		return personaOffsets[u] + i;
	};

	// Connect personas of each node
	double weightFactor = std::stod(parameters.at("personaEdgeWeightFactor"));
	G.forNodes([&](node u){
		for (auto edge : personaEdges[u]) {
			personaGraph.addEdge(getPersona(u, std::get<0>(edge)), getPersona(u, std::get<1>(edge)),
			                     std::get<2>(edge) * weightFactor);
		}
	});

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
	for (const auto& edges : personaEdges)
		internalPersonaEdges += edges.size();
	assert(personaGraph.numberOfEdges() == G.numberOfEdges() + internalPersonaEdges);
	// check that no isolated nodes were added
	auto numIsolatedNodes = [](const Graph &G) {
		ConnectedComponents compsAlgo(G);
		compsAlgo.run();
		auto comps = compsAlgo.getComponentSizes();
		count isolated = 0;
		for (const auto &x  : comps) {
			if (x.second == 1)
				++isolated;
			assert(x.second != 0);
			assert(x.second != 1);
		}
		return isolated;
	};
	auto gEdges = G.edges();
	auto pEdges = personaGraph.edges();
	//auto iso1 = numIsolatedNodes(G);
//	auto iso2 = numIsolatedNodes(personaGraph);
//	assert(iso2 == 0);
#endif
}

void EgoSplitting::createPersonaClustering() {
	personaPartition = globalClusterAlgo(personaGraph);
	assert(personaPartition.upperBound() <= personaGraph.upperNodeIdBound());
}

void EgoSplitting::createCover() {
	// Create cover from persona partition
	cover = Cover(G.upperNodeIdBound());
	personaPartition.compact();
	cover.setUpperBound(personaPartition.upperBound());
	G.forNodes([&](node u) {
		for (index i = personaOffsets[u]; i < personaOffsets[u + 1]; ++i) {
			cover.addToSubset(personaPartition.subsetOf(i), u);
		}
	});

	// Discard communities of size 4 or less
	count min_size = 5;
	std::vector<std::vector<node>> communities(cover.upperBound());
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

Partition
EgoSplitting::createGroundTruthPartition(Graph &egoGraph, NodeMapping &mapping, node egoNode) const {
	auto truthComms = groundTruth.subsetsOf(egoNode);
	index subset = groundTruth.upperBound();
	Partition part(subset + egoGraph.upperNodeIdBound());
	egoGraph.forNodes([&](node v){
		auto comms = groundTruth.subsetsOf(mapping.global(v));
		std::vector<node> overlap(truthComms.size());
		std::set_intersection(truthComms.begin(), truthComms.end(), comms.begin(), comms.end(),
		                      overlap.begin());
		if (!overlap.empty()) {
			part.addToSubset(overlap[0], v);
		} else {
//			part.addToSubset(subset, v);
//			++subset;
		}
	});
	return part;
}

Cover EgoSplitting::getCover() {
	return cover;
}

std::map<std::string, double> EgoSplitting::getTimings() {
	return timings;
}

std::map<std::string, double> EgoSplitting::getExecutionInfo() {
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
