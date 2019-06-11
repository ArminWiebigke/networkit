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
#include "../oslom/Stochastics.h"

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
	directedEdges = AdjacencyArray(G);
	Stochastics::init(G.numberOfEdges());

	parameters["storeEgoNet"] = "No";
	parameters["addEgoNode"] = "No";
	parameters["partitionFromGroundTruth"] = "No";
	parameters["extendFromPartitionIterations"] = "2";
	parameters["extendStrategySecond"] = "significance";

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
	parameters["processEgoNet"] = "extend";
	parameters["extendRandom"] = "No";
	parameters["extendOverDirected"] = "No";
	parameters["addNodesFactor"] = "4";
	parameters["addNodesExponent"] = "0.6";
	parameters["edgesBetweenNeigNeig"] = "Yes";
	parameters["keepOnlyTriangles"] = "No";
	parameters["minNodeDegree"] = "0";
	parameters["triangleThreshold"] = "0";

	// Parameters for edgeScores
	parameters["extendStrategy"] = "edgeScore";
	parameters["scoreStrategy"] = "score";

	// Test parameters
//	parameters["storeEgoNet"] = "Yes";
//	parameters["connectPersonas"] = "No";
//	parameters["processEgoNet"] = "extend";
//	parameters["extendStrategy"] = "significance";
//	parameters["scoreStrategy"] = "score";
	parameters["maxSignificance"] = "0.1";

}

void EgoSplitting::run() {
	Aux::SignalHandler handler;
	timings.clear();
	Aux::Timer timer;
	timer.start();

	if (std::stod("0.3") == 0)
		throw std::runtime_error("Can't convert numbers because of wrong locale!");
	G.forNodes([&](node u) {
		if (G.hasEdge(u, u))
			throw std::runtime_error("No self-loops allowed!");
	});


//	edgeScoreGraph = weightedEdgesGraph(G);

	INFO("create EgoNets");
	createEgoNets();
	addTime(timer, "1    create EgoNets");
	handler.assureRunning();

	INFO("split into Personas");
	splitIntoPersonas();
	addTime(timer, "2    split Personas");
	handler.assureRunning();

	INFO("connect Personas");
	connectPersonas();
	addTime(timer, "3    connect Personas");
	handler.assureRunning();

	INFO("create Persona Clustering");
	timer.start();
	createPersonaClustering();
	timer.stop();
	addTime(timer, "4    Persona Clustering");
	handler.assureRunning();

	INFO("create Cover");
	timer.start();
	createCover();
	timer.stop();
	addTime(timer, "5    create Cover");

	hasRun = true;
}

std::string EgoSplitting::toString() const {
	return "EgoSplitting";
}

//Graph EgoSplitting::extendEgoNet(node u, Graph egoGraph, const NodeMapping &nodeMapping) {
//	double addNodesFactor = std::stod(parameters.at("addNodesFactor"));
//	double addNodesExponent = std::stod(parameters.at("addNodesExponent"));
//	count extendNodeCnt = std::ceil(
//			addNodesFactor * std::pow(egoGraph.numberOfNodes(), addNodesExponent));
//	if (parameters.at("processEgoNet") == "none") {
//		// nothing to do
//	} else if (parameters.at("processEgoNet") == "extend") {
//		extendEgoNet(egoGraph, u, nodeMapping, extendNodeCnt);
//	} else {
//		throw std::runtime_error("Missing strategy to extend ego-net");
//	}
//	// TODO: score basierend auf edges / degree des Kandidaten
//	// TODO: Überlappende Partitionierung
//	return egoGraph;
//}

Partition EgoSplitting::partitionEgoNet(node u, const Graph &egoGraph,
                                        const NodeMapping &nodeMapping) const {
	// Partition
	Partition egoPartition;
	if (parameters.at("partitionFromGroundTruth") == "Yes")
		egoPartition = createGroundTruthPartition(egoGraph, nodeMapping, u);
	else if (egoGraph.numberOfEdges() > 0) {
		egoPartition = localClusterAlgo(egoGraph);
	} else {
		egoPartition = Partition(egoGraph.upperNodeIdBound());
		egoPartition.allToSingletons();
	}
	egoPartition.compact();
	return egoPartition;
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
		INFO("Add neighbors");
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
		addTime(timer, "11    Find nodes");


		/******************************************************************************************
		 **                             Triangle Search for Edges                                **
		 ******************************************************************************************/
		INFO("Add edges");
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
		addTime(timer, "12    Neighbor Triangle Search");

		/*
		 * TODO: Zwei mal partitionieren:
		 * Beim ersten Partitionieren: Extende maximal (niedrieger Wert?) viele Knoten, einfach nur
		 * nach Score sortiert.
		 * Beim zweiten Partitionieren: Extende nur Knoten, die zu einer Partition oder zu einer
		 * Vereinigung von Partitionen signifikant sind. (oder maximal viele)
		 *
		 * TODO: Vereinigung von Partitionen überprüfen:
		 * Verbinde nur Partitionen die auch verbunden sind.
		 * Sortiere nach Anzahl Kanten. Überprüfe Partitionen wie folgt, Abbruch falls ein Kandidat
		 * sig. ist.
		 * Berechne Sig. für erste Part.
		 * Berechne Sig. für nächste Part., füge sie in sortierte Liste der Sig.en an Stelle x ein.
		 * Berechne Sig. für Vereinigung der x besten Part.
		 *
		 *
		 */

		/******************************************************************************************
		 **                          Extend and Partition EgoNet                                 **
		 ******************************************************************************************/
		INFO("Extend EgoNet");
		Graph egoGraphBase(egoGraph);
		NodeMapping nodeMappingBase(nodeMapping);
		Partition egoPartition;
		count extendIterations = std::stoi(parameters.at("extendFromPartitionIterations"));
		for (count i = 0; i < extendIterations; ++i) {
			if (i > 0) {
				egoGraph = egoGraphBase;
				nodeMapping = nodeMappingBase;
			}

			extendEgoNet(u, egoGraph, nodeMapping, egoPartition);
			addTime(timer, "15)    Extend EgoNet");

			egoPartition = partitionEgoNet(u, egoGraph, nodeMapping);
			addTime(timer, "16)    Cluster EgoNet");
		}


		/******************************************************************************************
		 **                                 Connect Personas                                     **
		 ******************************************************************************************/
		INFO("Connect Personas of one node");
		if (parameters.at("connectPersonas") == "Yes")
			personaEdges[u] = connectEgoPartitionPersonas(egoGraph, egoPartition);
		addTime(timer, "17)    Connect Personas");



		/******************************************************************************************
		 **                                 Store EgoNet                                         **
		 ******************************************************************************************/
		/* Store EgoNet for evaluation with global node IDs  */
		if (parameters.at("storeEgoNet") == "Yes") {
			INFO("Store EgoNet");
			// Get EgoNet with gloabl node ids
			Graph egoNetGraph = Graph(G.upperNodeIdBound(), egoGraph.isWeighted());
			egoGraph.forEdges([&](node v, node w, edgeweight weight) {
				egoNetGraph.addEdge(nodeMapping.global(v), nodeMapping.global(w), weight);
			});
			for (node v : egoNetGraph.nodes()) {
				if (!egoGraph.hasNode(nodeMapping.local(v)))
					egoNetGraph.removeNode(v);
			}
			if (egoGraph.isWeighted() && egoNetGraph.numberOfNodes() >= 2)
				egoNetGraph.addEdge(egoNetGraph.nodes()[0], egoNetGraph.nodes()[1],
				                    0.000001);
			egoNets[u] = egoNetGraph;
		}
		addTime(timer, "1b)    Copy EgoNet");


		/******************************************************************************************
		 **                                 Build EgoNet Map                                     **
		 ******************************************************************************************/
		INFO("Build EgoNet Partition Map");
		// Insert nodes into ego-net data structure
		for (node i : egoGraph.nodes()) {
			egoNetPartitions[u].emplace(nodeMapping.global(i), egoPartition.subsetOf(i));
		}
		assert(egoNetPartitions[u].size() == egoGraph.numberOfNodes());
		egoNetPartitionCounts[u] = egoPartition.numberOfSubsets();
		addTime(timer, "1c)    EgoNet subsets");


		/******************************************************************************************
		 **                                     Cleanup                                          **
		 ******************************************************************************************/
		nodeMapping.reset();
		addTime(timer, "1x)    Clean up");
	});
}

std::vector<std::tuple<index, index, edgeweight>>
EgoSplitting::connectEgoPartitionPersonas(const Graph &egoGraph,
                                          const Partition &egoPartition) const {
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

	std::string strategy = parameters.at("connectPersonasStrat");
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
		throw std::runtime_error(strategy + " is not a valid strategy to connect personas!");
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
EgoSplitting::extendEgoNet(node u, Graph &egoGraph, NodeMapping &nodeMapping,
                           Partition &basePartition) const {
	if (parameters.at("processEgoNet") == "none") {
		// nothing to do
		return;
	} else if (parameters.at("processEgoNet") == "extend") {
		// continue below
	} else {
		throw std::runtime_error("Missing strategy to extend ego-net");
	}

	Aux::Timer timer;
	timer.start();
	double addNodesFactor = std::stod(parameters.at("addNodesFactor"));
	double addNodesExponent = std::stod(parameters.at("addNodesExponent"));
	count extendNodeCnt = std::ceil(
			addNodesFactor * std::pow(egoGraph.numberOfNodes(), addNodesExponent));
	const count directNeighborsCnt = nodeMapping.nodeCount();
	std::vector<std::pair<node, double>> nodeScores; // node and its score
	bool useBasePartition = basePartition.numberOfSubsets() > 0;
	if (!useBasePartition) {
		basePartition = Partition(egoGraph.upperNodeIdBound());
		basePartition.allToOnePartition();
	}
	addTime(timer, "150    Setup");

	/**********************************************************************************************
	 **                           Get node candidates with scores                                **
	 **********************************************************************************************/
	std::string extendStrategy = parameters.at("extendStrategy");
	if (useBasePartition)
		extendStrategy = parameters.at("extendStrategySecond");

	std::vector<std::set<node>> triangleEdges;
	if (extendStrategy == "edgeScore")
		nodeScores = scoreEdgeCount(u, nodeMapping);
	else if (extendStrategy == "triangles")
		nodeScores = scoreTriangles(u, nodeMapping, triangleEdges, egoGraph);
	else if (extendStrategy == "significance")
		nodeScores = scoreSignificance(u, nodeMapping, egoGraph, basePartition);
	else
		throw std::runtime_error(extendStrategy + " is not a valid strategy to extend Ego-Net!");

#ifndef NDEBUG
	std::set<node> candidates;
	for (auto pair : nodeScores) {
		node v = pair.first;
		assert(candidates.count(v) == 0);
		candidates.insert(v);
	}
	assert(candidates.size() == nodeScores.size());
#endif

	// Remove nodes with score zero
	for (size_t i = 0; i < nodeScores.size(); ++i) {
		if (nodeScores[i].second <= 0.0) {
			nodeScores[i] = nodeScores.back();
			nodeScores.pop_back();
			--i; // Check node at current position again
		}
	}
	addTime(timer, "151 Get candidates " + extendStrategy);
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
		nodeMapping.addNode(v);
		egoGraph.addNode();
	}
	Graph onlyNeighborsGraph(egoGraph);
	addTime(timer, "153 Add nodes");


	/**********************************************************************************************
	 **                                  Add edges to result                                     **
	 **********************************************************************************************/
	for (node v : nodeMapping.globalNodes()) {
		directedEdges.forEdgesOf(v, [&](node, node w, edgeweight weight) {
			if (nodeMapping.isMapped(w)) {
				node v_loc = nodeMapping.local(v);
				node w_loc = nodeMapping.local(w);
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
	addTime(timer, "155 Add edges");


	/**********************************************************************************************
	 **                         Use only triangle edges (optional)                               **
	 **********************************************************************************************/
	if (parameters.at("keepOnlyTriangles") == "Yes") {
		AdjacencyArray directedEgoGraph(egoGraph);
		auto foundTriangleWork = [&](node v, node w, node x) {
			node v_loc = nodeMapping.local(v);
			node w_loc = nodeMapping.local(w);
			node x_loc = nodeMapping.local(x);
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

	count minDegree = stoi(parameters.at("minNodeDegree"));
	removeLowDegreeNodes(egoGraph, minDegree, directNeighborsCnt);
	removeLowTriangleCntNodes(egoGraph,
	                          directNeighborsCnt); // TODO: Bevorzuge high degree nodes (?)
	removeLowDegreeNodes(egoGraph, minDegree, directNeighborsCnt);
	addTime(timer, "157 Remove low degree nodes");
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
			if (v >= directNeighborsCnt &&
			    triangleCount[v] * 1.0 / egoGraph.degree(v) < triangleThreshold) {
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
	throw std::runtime_error(scoreStrategy + " is not a valid score strategy!");
}

std::vector<std::pair<node, double>>
EgoSplitting::scoreEdgeCount(node u, const NodeMapping &neighbors) const {
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
                             Graph const &egoGraph) const {
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

void EgoSplitting::addTime(Aux::Timer &timer, const std::string &name) const {
	timer.stop();
	double elapsed = timer.elapsedNanoseconds();
	timings[name] += elapsed;
	timer.start();
}


std::vector<std::pair<node, double>>
EgoSplitting::scoreSignificance(node u, const NodeMapping &egoMapping, Graph const &egoGraph,
                                const Partition &basePartition) const {
	Aux::Timer timer;
	timer.start();
	ParallelPartitionCoarsening coarsening(egoGraph, basePartition);
	coarsening.run();
	Graph coarseGraph = coarsening.getCoarseGraph();
	auto coarseToEgo = coarsening.getCoarseToFineNodeMapping();
	auto egoToCoarse = coarsening.getFineToCoarseNodeMapping();

	std::vector<count> coarseSizes(coarseGraph.upperNodeIdBound());
	coarseGraph.forNodes([&](node p) {
		coarseSizes[p] = coarseToEgo[p].size();
	});

	addTime(timer, "151s0    Coarsening");
	// TODO: "Remove" egoNode from graph for calculation

	// Insert a node that represents all external nodes (that are not candidates)
	node externalNode = coarseGraph.addNode();
	coarseSizes.push_back(G.numberOfNodes() - egoGraph.numberOfNodes());

	// NodeMapping zwischen global und coarse für die Kandidaten
	NodeMapping coarseMapping(G); // Maps global graph to coarseGraph
	coarseGraph.forNodes([&](node p) {
		coarseMapping.addDummy();
	});
	assert(coarseMapping.nodeCount() == coarseGraph.upperNodeIdBound());

	// Get candidates
	std::vector<std::map<node, double>> edgeScores(G.upperNodeIdBound());
	std::set<node> candidates;
	std::vector<node> directNeighbors = egoMapping.globalNodes();
	auto isDirectNeighbor = [&](node x) {
		return egoMapping.isMapped(x);
	};
	auto countEdges = [&](node v, node w, edgeweight weight) {
		if (!isDirectNeighbor(w)) {
			assert(!coarseMapping.isMapped(w));
			if (edgeScores[w].empty())
				candidates.insert(w);
			edgeScores[w][egoToCoarse[egoMapping.local(v)]] += weight;
		}
	};
	for (node v : directNeighbors) {
		if (parameters.at("extendOverDirected") == "Yes")
			directedEdges.forEdgesOf(v, countEdges);
		else
			G.forEdgesOf(v, countEdges);
	}
	candidates.erase(u); // Remove egoNode as candidate

	// Füge Kandidaten in Graph ein
	for (node v : candidates) {
		assert(!coarseMapping.isMapped(v));
		coarseMapping.addNode(v);
		coarseGraph.addNode();
		for (auto p : edgeScores[v]) {
			assert(p.second <= coarseSizes[p.first]);
			coarseGraph.addEdge(coarseMapping.local(v), p.first, p.second);
		}
	}
	addTime(timer, "151s2    Add candidates");


	// Check
	auto check = [&]() {
		for (node v : candidates) {
			std::map<node, count> groupEdges;
			G.forNeighborsOf(v, [&](node w) {
				if (egoMapping.isMapped(v))
					groupEdges[egoToCoarse[w]] += 1;
			});
			for (auto pair : groupEdges)
				assert(pair.second == coarseGraph.weight(v, pair.first));
		}
		assert(coarseGraph.upperNodeIdBound() == coarseMapping.nodeCount());
		return true;
	};
	assert(check());

	edgeweight totalWeight = coarseGraph.totalEdgeWeight();
	count numEgoEdges = G.numberOfEdges();
	count missingEdges = numEgoEdges - totalWeight;
	coarseGraph.addEdge(externalNode, externalNode, missingEdges);

	return calcSignficance(externalNode, coarseGraph, coarseMapping, coarseSizes);
}

std::vector<std::pair<node, double>>
EgoSplitting::calcSignficance(node externalNode, const Graph &coarseGraph,
                              const NodeMapping &coarseMapping,
                              const std::vector<count> &coarseSizes) const {
	Aux::Timer timer;
	timer.start();
	/* Calculate score for each candidate */
	// Sort candidates by number of edges
	std::vector<std::pair<count, node>> candidatesSorted;
//	nodes_and_edges.reserve(candidates.size());
	for (node v = externalNode + 1; v < coarseGraph.upperNodeIdBound(); ++v) {
		count edgesIntoEgo = coarseGraph.weightedDegree(v);
		candidatesSorted.emplace_back(edgesIntoEgo, coarseMapping.global(v));
	}
	std::sort(candidatesSorted.rbegin(), candidatesSorted.rend());
	addTime(timer, "151s6    Sort candidates");

	// Discard nodes with less than 3 edges
	for (count i = 0; i < candidatesSorted.size(); ++i) {
		if (candidatesSorted[i].first < 3) {
			candidatesSorted.resize(i);
			break;
		}
	}
	addTime(timer, "151s8    Discard candidates");

	std::vector<count> group_total(externalNode);
	std::vector<count> group_outgoing(externalNode);
	std::vector<count> external_stubs(externalNode);
	std::vector<count> external_nodes(externalNode);
	for (count p = 0; p < externalNode; ++p) {
		if (!coarseGraph.hasNode(p))
			continue;

		group_total[p] = (int) coarseGraph.weightedDegree(p);
		group_outgoing[p] = coarseGraph.weightedDegree(p) - 2 * (int) coarseGraph.weight(p, p);
		external_stubs[p] = G.numberOfEdges() * 2 - group_total[p];
		external_nodes[p] = G.numberOfNodes() - coarseSizes[p];
	}

	std::vector<std::pair<node, double>> nodeScores;
	double max_significance = std::stod(parameters.at("maxSignificance"));
	for (auto pair : candidatesSorted) {
		node v = pair.second;
		addTime(timer, "151sg    Loop");
		count node_degree = G.degree(v);
		node loc_v = coarseMapping.local(v);
		addTime(timer, "151s9    Get degree and local");
		// TODO: Check group merge
		std::vector<std::pair<double, node>> groupEdges;
		auto groups = coarseGraph.neighbors(loc_v);
		for (node p : groups) {
			if (p >= externalNode)
				continue;
			edgeweight numEdges = coarseGraph.weight(loc_v, p);
			groupEdges.emplace_back(numEdges, p);
			assert(numEdges <= coarseSizes[p]);
		}
		std::sort(groupEdges.rbegin(), groupEdges.rend());

		auto addIfSignificant = [&](double significance) {
			if (significance <= max_significance) {
				double score = 1 - significance;
				score = normalizeScore(v, score);
				nodeScores.emplace_back(v, score);
				return true;
			}
			return false;
		};
		bool added = false;
		count counter = 0;
		for (auto it : groupEdges) {
			if (++counter > 3)
				break;
			addTime(timer, "151sg    Loop");
			count numEdges = (int) it.first;
			node p = it.second;
			addTime(timer, "151sa    Get num edges");
			if (numEdges < 3)
				break;
			double significance =
					Stochastics::calc_score(node_degree, (int) numEdges, group_outgoing[p],
					                        external_stubs[p], external_nodes[p]);
			added = addIfSignificant(significance);
			addTime(timer, "151sc    Calc significance");
			if (added)
				break;
		}
		if (added)
			continue;

		// TODO: Nimm nicht den besten sondern den x-besten (z.B. 10% von EgoNet) für die
		// Ordered Statistics. Besser als maxSignificance zu skalieren.

		auto it = groupEdges.begin();
		node bestGroup = it->second;
		std::set<node> mergedGroups{bestGroup};
		count group_total_stubs = group_total[bestGroup];
		count group_out_stubs = group_outgoing[bestGroup];
		count ext_stubs = external_stubs[bestGroup];
		count ext_nodes = external_nodes[bestGroup];
		count numEdges = (int) it->first;

		auto check = [&]() {
			// Check
			count ext_nodes_check = G.numberOfNodes();
			for (node w : mergedGroups)
				ext_nodes_check -= coarseSizes[w];
			assert(ext_nodes == ext_nodes_check);
			count group_out_stubs_check = 0;
			count group_total_stubs_check = 0;
			for (node w : mergedGroups) {
				coarseGraph.forEdgesOf(w, [&](node, node ext, edgeweight weight) {
					group_total_stubs_check += (int) weight;
					if (mergedGroups.count(ext) == 0)
						group_out_stubs_check += (int) weight;
					if (w == ext)
						group_total_stubs_check += (int) weight;
				});
			}
			assert(group_total_stubs_check == group_total_stubs);
			assert(group_out_stubs == group_out_stubs_check);
			count ext_stubs_check = 2 * G.numberOfEdges() - group_total_stubs_check;
			assert(ext_stubs_check == ext_stubs);

			assert(ext_stubs_check >= group_out_stubs_check);


			count sum = 0;
			coarseGraph.forEdges([&](node n1, node n2, edgeweight weight) {
				sum += 2 * (int) weight;
			});
			count g_in = 0;
			count g_out = 0;
			count ext = 0;
			coarseGraph.forEdges([&](node n1, node n2, edgeweight weight) {
				count id = mergedGroups.count(n1) + mergedGroups.count(n2);
				if (id == 0) {
					ext += 2 * (int) weight;
				} else if (id == 1) {
					ext += (int) weight;
					g_out += (int) weight;
				} else if (id == 2) {
					g_in += 2 * (int) weight;
				}
			});
			assert(sum == g_in + g_out + ext);
			assert(g_out <= ext);
			assert(g_out == group_out_stubs);
			assert(ext == ext_stubs);
			return true;
		};

		assert(check());
		if (parameters.at("sortGroupsDensity") == "Yes") {
			std::sort(groupEdges.begin(), groupEdges.end(), [&](std::pair<double, node> a,
			                                                    std::pair<double, node> b) {
				return (a.first / coarseSizes[a.second]) > b.first / coarseSizes[b.second];
			});
		}

		for (++it; it < groupEdges.end(); ++it) {
			node group = it->second;
			// Merge group
			ext_stubs -= group_total[group];
			group_total_stubs += group_total[group];
			group_out_stubs += group_outgoing[group];
			for (node w : mergedGroups) {
				group_out_stubs -= 2 * (int) coarseGraph.weight(w, group); // internal stubs
			}
			ext_nodes -= coarseSizes[group];
			numEdges += (int) it->first;
			mergedGroups.insert(group);

			assert(check());

			// Calculate new significance
			double significance =
					Stochastics::calc_score(node_degree, numEdges, group_out_stubs,
					                        ext_stubs, ext_nodes);
			added = addIfSignificant(significance);
			if (added) {
				std::cout << mergedGroups.size() << " merged groups " << nodeScores.back().second
				          << " (" << numEdges << "/" << node_degree << ") -> "
				          << G.numberOfNodes() - ext_nodes
				          << std::endl;
				break;
			}
			if (mergedGroups.size() == 4)
				break;
		}
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
	G.forNodes([&](node u) {
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
	for (const auto &edges : personaEdges)
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
EgoSplitting::createGroundTruthPartition(const Graph &egoGraph, const NodeMapping &mapping,
                                         node egoNode) const {
	auto truthComms = groundTruth.subsetsOf(egoNode);
	index subset = groundTruth.upperBound();
	Partition part(subset + egoGraph.upperNodeIdBound());
	egoGraph.forNodes([&](node v) {
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
