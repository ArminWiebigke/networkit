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
#include "EgoSplitting.h"
#include "PLP.h"
#include "CoverF1Similarity.h"

namespace NetworKit {

EgoSplitting::EgoSplitting(const Graph &G)
		: G(G) {
	std::function<Partition(Graph &)> clusterAlgo = [](Graph &G) {
		PLP algo(G, 1, 20);
		algo.run();
		return algo.getPartition();
	};
	localClusterAlgo = clusterAlgo;
	globalClusterAlgo = clusterAlgo;
	init();
}

EgoSplitting::EgoSplitting(const Graph &G,
                           std::function<Partition(Graph &)> clusterAlgo)
		: EgoSplitting(G, clusterAlgo, clusterAlgo) {
}

EgoSplitting::EgoSplitting(const Graph &G,
                           std::function<Partition(Graph &)> localClusterAlgo,
                           std::function<Partition(Graph &)> globalClusterAlgo)
		: G(G),
		  localClusterAlgo(std::move(localClusterAlgo)),
		  globalClusterAlgo(std::move(globalClusterAlgo)) {
	init();
}

void EgoSplitting::init() {
	egoNets.resize(G.upperNodeIdBound());
	egoNetPartitions.resize(G.upperNodeIdBound());
	egoNetPartitionCounts.resize(G.upperNodeIdBound(), 0);
	personaOffsets.resize(G.upperNodeIdBound() + 1, 0);
	directedEdges = AdjacencyArray(G);

	parameters["storeEgoNet"] = "No";
	parameters["addEgoNode"] = "No";

	// Parameters for weighted edges
	parameters["weightFactor"] = "1";
	parameters["weightOffset"] = "0";

	// Parameters for ego-net extension
	parameters["processEgoNet"] = "extend";
	parameters["addNodesFactor"] = "1";
	parameters["addNodesExponent"] = "0.6";

	// Parameters for edgeScores
	parameters["scoreStrategy"] = "count";
	parameters["extendRandom"] = "No";
	parameters["keepOnlyTriangles"] = "No";


	// Test parameters
	parameters["processEgoNet"] = "extend";
	parameters["extendRandom"] = "No";
	parameters["extendStrategy"] = "triangles";
	parameters["scoreStrategy"] = "score";
}

void EgoSplitting::run() {
	Aux::SignalHandler handler;
	Aux::Timer timer;
	timings.clear();

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
	double weightFactor = std::stod(parameters["weightFactor"]);
	double weightOffset = std::stod(parameters["weightOffset"]);
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
		DEBUG("Create EgoNet for Node ", u, "/", G.upperNodeIdBound());

		/******************************************************************************************
		 **                                Add Neighbor Nodes                                    **
		 ******************************************************************************************/
		timer.start();
		count degree = G.degree(u);
		Graph egoGraph(degree, true);

		if (parameters["addEgoNode"] == "Yes") {
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
			if (parameters["addEgoNode"] == "Yes") {
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
		double addNodesFactor = std::stod(parameters["addNodesFactor"]);
		double addNodesExponent = std::stod(parameters["addNodesExponent"]);
		count extendNodeCnt = std::ceil(addNodesFactor * std::pow(degree, addNodesExponent));
		if (parameters["processEgoNet"] == "none") {
			// nothing to do
		} else if (parameters["processEgoNet"] == "extend") {
			extendEgoNet(egoGraph, u, nodeMapping, extendNodeCnt);
		} else if (parameters["processEgoNet"] == "consensus") {
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
		if (egoGraph.numberOfEdges() > 0) {
			Graph egoCopy(egoGraph, false, egoGraph.isDirected());
			egoPartition = localClusterAlgo(egoCopy); // Python moves the graph so we need a copy
		} else {
			egoPartition = Partition(degree);
			egoPartition.allToSingletons();
		}
		egoPartition.compact();
		timer.stop();
		timings["1d)    Cluster EgoNet"] += timer.elapsedMicroseconds();


		/******************************************************************************************
		 **                                  Store EgoNet                                        **
		 ******************************************************************************************/
		timer.start();
		if (parameters["storeEgoNet"] == "Yes") {
			// Get EgoNet with gloabl node ids
			egoNets[u] = Graph(G.upperNodeIdBound(), true);
			egoGraph.forEdges([&](node v, node w, edgeweight weight) {
				egoNets[u].addEdge(nodeMapping.global(v), nodeMapping.global(w), weight);
			});
			for (node v : egoNets[u].nodes()) {
				if (!egoGraph.hasNode(nodeMapping.local(v)))
					egoNets[u].removeNode(v);
			}
			egoNets[u].addEdge(egoNets[u].nodes()[0], egoNets[u].nodes()[1], 0.000001);
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

void EgoSplitting::consensusWeighting(Graph &egoGraph, node u, const NodeMapping &nodeMapping,
                                      count extendNodeCnt) {
	if (parameters["extendRandom"] != "Yes")
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

void
EgoSplitting::extendEgoNet(Graph &egoGraph, node u, NodeMapping &neighbors, count extendNodeCnt) {
	count directNeighborsCnt = neighbors.nodeCount();
	std::vector<std::pair<node, double>> nodeScores; // node and its score

	/**********************************************************************************************
	 **                           Get node candidates with scores                                **
	 **********************************************************************************************/
	const auto extendStrategy = parameters["extendStrategy"];
	std::vector<std::set<node>> triangleEdges;
	if (extendStrategy == "edgeScore")
		nodeScores = scoreEdgeCount(u, neighbors);
	else if (extendStrategy == "triangles")
		nodeScores = scoreTriangles(u, neighbors, triangleEdges);
	else
		throw std::runtime_error("No valid strategy to extend Ego-Net");

	// Remove nodes with score zero
	for (size_t i = 0; i < nodeScores.size(); ++i) {
		if (nodeScores[i].second <= 0.0) {
			nodeScores[i] = nodeScores.back();
			nodeScores.pop_back();
		}
	}
	// TODO: Use only nodes with score > 0.1 (or something), use all nodes with score > 0.5
	// TODO: Score need to be normalized to (0, 1) for this to work


	/**********************************************************************************************
	 **                                  Add nodes to result                                     **
	 **********************************************************************************************/
	std::vector<node> nodesToAdd;
	if (parameters["extendRandom"] == "Yes") {
		// Take nodes randomly, the chance to pick a node is proportional to its score
		std::vector<double> random_weights;
		for (auto pair : nodeScores)
			random_weights.push_back(pair.second);

		for (index i = 0; i < extendNodeCnt && !nodeScores.empty(); ++i) {
			std::discrete_distribution<index> dist(random_weights.begin(), random_weights.end());
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
	for (auto v : nodesToAdd) {
		neighbors.addNode(v);
		egoGraph.addNode();
	}


	/**********************************************************************************************
	 **                                  Add edges to result                                     **
	 **********************************************************************************************/
	if (parameters["keepOnlyTriangles"] == "Yes") {
		for (node v : nodesToAdd) {
			for (node w : triangleEdges[v])
				egoGraph.addEdge(neighbors.local(v), neighbors.local(w), 1);
		}
	} else {
		for (node v : neighbors.globalNodes()) {
			G.forEdgesOf(v, [&](node, node w, edgeweight weight) {
				if (neighbors.isMapped(w)) {
					// Discard edges between neighbors of neighbors
					if (parameters["edgesBetweenNeigNeig"] == "No"
					    && neighbors.local(v) >= directNeighborsCnt
					    && neighbors.local(w) >= directNeighborsCnt) {
						return;
					}
					if (!egoGraph.hasEdge(neighbors.local(v), neighbors.local(w)))
						egoGraph.addEdge(neighbors.local(v), neighbors.local(w), weight);
				}
			});
		}
	}
}

double EgoSplitting::normalizeScore(node v, double score) {
	std::string scoreStrategy = parameters["scoreStrategy"];
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
	for (node v : directNeighbors) {
		G.forEdgesOf(v, [&](node, node w, edgeweight weight) {
			if (!isDirectNeighbor(w) &&
			    w != u) { // TODO: Adding the node u to the ego-net seems to improve the result
				if (edgeScores[w].empty())
					secondNeighbors.push_back(w);
				edgeScores[w].push_back(weight);
			}
		});
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
                             std::vector<std::set<node>> &triangleEdges) {
	std::vector<int> numTriangles(G.upperNodeIdBound(), -1);
	std::vector<node> secondNeighbors;
	NodeMapping directAndSecondNeighbors = neighbors;
	std::vector<node> directNeighbors = neighbors.globalNodes();
	auto isDirectNeighbor = [&](node x) {
		return neighbors.isMapped(x);
	};
	for (node v : directNeighbors) {
		G.forEdgesOf(v, [&](node, node w, edgeweight weight) {
			if (!isDirectNeighbor(w) && w != u) {
				if (numTriangles[w] == -1) {
					secondNeighbors.push_back(w);
					directAndSecondNeighbors.addNode(w);
				}
				numTriangles[w] = 0;
			}
		});
	}
	auto getSingleSecondNeighbor = [&](node node1, node node2, node node3){
		int b1 = numTriangles[node1] != -1;
		int b2 = numTriangles[node2] != -1;
		int b3 = numTriangles[node3] != -1;
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
	triangleEdges.resize(G.upperNodeIdBound());
	auto countTriangles = [&](node v){
		std::vector<int> currentNeighbors(G.upperNodeIdBound(), 0);
		directedEdges.forEdgesOf(v, [&](node, node w, edgeweight weight1) {
			currentNeighbors[w] = 1;
		});
		directedEdges.forEdgesOf(v, [&](node, node w, edgeweight weight1) {
			if (!directAndSecondNeighbors.isMapped(w))
				return;
			directedEdges.forEdgesOf(w, [&](node, node x, edgeweight weight2) {
				if (!directAndSecondNeighbors.isMapped(x) || currentNeighbors[x] == 0)
					return;
				// We found a triangle v-w-x
				assert(G.hasEdge(v, w) && G.hasEdge(w, x) && G.hasEdge(x, v));
				// Only count triangles if exactly one node is a second neighbor
				node secondNeighbor = getSingleSecondNeighbor(v, w, x);
				if (secondNeighbor != none) {
					assert(!G.hasEdge(u, secondNeighbor));
					++numTriangles[secondNeighbor];
					for (node neighborNode : {v, w, x}){
						if (secondNeighbor != neighborNode) {
							assert(G.hasEdge(u, neighborNode));
							triangleEdges[secondNeighbor].insert(neighborNode);
						}
					}
				}
			});
		});
	};
	for (node v : directNeighbors) {
		countTriangles(v);
	}
	for (node v : secondNeighbors) {
		countTriangles(v);
	}

	std::vector<std::pair<node, double>> nodeScores;
	nodeScores.reserve(secondNeighbors.size());
	for (node v : secondNeighbors) {
		double score = numTriangles[v]; // TODO: Maybe use number of triangle edges as score?
		// TODO: num_triangles / num_triangle_edges : Sinnvolle Metrik?
		double normalizedScore = normalizeScore(v, score);
		nodeScores.emplace_back(v, normalizedScore);
	}
	return nodeScores;
}


//	Graph extendedGraph = Graph(egoGraph.upperNodeIdBound());
//	for (node v : extendedGraph.nodes()) {
//		if (!egoGraph.hasNode(v))
//			extendedGraph.removeNode(v);
//	}
//
//	AdjacencyArray directedEgoGraph(egoGraph);
//	std::vector<double> neighborEdgeWeight(egoGraph.upperNodeIdBound(), 0.0);
//	egoGraph.forNodes([&](node v) {
//		// Mark neighbors
//		directedEgoGraph.forEdgesOf(v, [&](node, node w, edgeweight weight) {
//			neighborEdgeWeight[w] = weight;
//		});
//
//		directedEgoGraph.forEdgesOf(v, [&](node, node w, edgeweight weight) {
//			// Always add edges between direct neighbors of u
//			if (v < directNeighborsCnt && w < directNeighborsCnt)
//				extendedGraph.addEdge(v, w, weight);
//			directedEgoGraph.forEdgesOf(w, [&](node, node x, edgeweight weight2) {
//				if (neighborEdgeWeight[x] != 0.0) {
//					// we have found a triangle v-w-x
//					auto tryAddEdge = [&](node u1, node u2, edgeweight edgeWeight) {
//						if ((u1 >= directNeighborsCnt || u2 >= directNeighborsCnt) &&
//						    !extendedGraph.hasEdge(u1, u2))
//							extendedGraph.addEdge(u1, u2, edgeWeight);
//					};
//					tryAddEdge(v, w, weight);
//					tryAddEdge(w, x, weight2);
//					tryAddEdge(v, x, neighborEdgeWeight[x]);
//				}
//			});
//		});
//
//		// Reset marked neighbors
//		directedEgoGraph.forEdgesOf(v, [&](node, node w, edgeweight weight) {
//			neighborEdgeWeight[w] = 0.0;
//		});
//	});

void EgoSplitting::splitIntoPersonas() {
	count sum = 0;
	for (index i = 0; i < G.upperNodeIdBound(); ++i) {
		personaOffsets[i] = sum;
		sum += egoNetPartitionCounts[i];
	}
	personaOffsets[G.upperNodeIdBound()] = sum;
	personaGraph = Graph(sum);
}

void EgoSplitting::connectPersonas() {
	auto getPersona = [&](node u, index i) {
		return personaOffsets[u] + i;
	};

	G.forEdges([&](node u, node v, edgeweight weight) {
		auto idx_u = egoNetPartitions[u].find(v);
		auto idx_v = egoNetPartitions[v].find(u);
		assert(idx_u != egoNetPartitions[u].end() && idx_v != egoNetPartitions[v].end());
		personaGraph.addEdge(getPersona(u, idx_u->second), getPersona(v, idx_v->second), weight);
	});

//	egoNetPartitions.clear();

#ifndef NDEBUG
	assert(personaGraph.numberOfEdges() == G.numberOfEdges());
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

Cover EgoSplitting::getCover() {
	return cover;
}

std::map<std::string, double> EgoSplitting::getTimings() {
	return timings;
}

std::map<std::string, double> EgoSplitting::getExecutionInfo() {
	return executionInfo;
}

std::vector<std::unordered_map<node, index>> EgoSplitting::getEgoNetPartitions() {
	return egoNetPartitions;
}

Graph EgoSplitting::getEgoNet(node u) {
	return egoNets[u];
}

void EgoSplitting::setParameters(std::map<std::string, std::string> const &new_parameters) {
	for (auto &x : new_parameters) {
		this->parameters[x.first] = x.second;
	}
}

} /* namespace NetworKit */
