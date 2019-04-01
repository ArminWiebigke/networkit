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

	parameters["storeEgoNet"] = "Yes";
	parameters["addEgoNode"] = "Yes";

	// Parameters for weighted edges
	parameters["weightFactor"] = "1";
	parameters["weightOffset"] = "0";
	parameters["discardThreshold"] = "0.0";  // Discard edges with smaller score

	// Parameters for ego-net extension
	parameters["extendStrategy"] = "edgeScores";
	parameters["addNodesFactor"] = "1";
	parameters["addNodesExponent"] = "0.8";

	// Parameters for simpleNN
	parameters["searchNeigOfNeighInDirectedGraph"] = "No";
	parameters["discardNeigOfNeigEdgesAtFirst"] = "Yes";
	parameters["discardNonTriangle"] = "Yes";
	parameters["minDegreeCleaning"] = "1";
	parameters["removeLowDegreeNeig"] = "Yes";
	parameters["edgesBetweenNeigNeig"] = "No";

	// Parameters for edgeScores
	parameters["scoreStrategy"] = "count";
}

void EgoSplitting::run() {
	Aux::SignalHandler handler;
	Aux::Timer timer;
	timings.clear();

	edgeScoreGraph = weightedEdgesGraph(G);

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
	auto toWeight = EdgeScoreAsWeight(inputGraph, localSimScores, false, weightOffset, weightFactor);
	Graph weightedGraph = toWeight.calculate();
	return weightedGraph;
}

void EgoSplitting::createEgoNets() {
//	AdjacencyArray directedEdges(G); // store each undirected edge as one directed edge
	AdjacencyArray directedEdges(edgeScoreGraph); // store each undirected edge as one directed edge
	// Assign IDs to the neighbours
	NodeMapping nodeMapping(G);

	Aux::SignalHandler handler;
	Aux::Timer timer;

	G.forNodes([&](node u) {
		handler.assureRunning();
		DEBUG("Create EgoNet for Node ", u, "/", G.upperNodeIdBound());
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

		// Extend ego-net
		double addNodesFactor = std::stod(parameters["addNodesFactor"]);
		double addNodesExponent = std::stod(parameters["addNodesExponent"]);
		count extendNodeCnt = std::ceil(addNodesFactor * std::pow(degree, addNodesExponent));
		if (parameters["extendStrategy"] == "none") {
			// nothing to do
		} else if (parameters["extendStrategy"] == "simpleNN") {
			extend_simpleNN(egoGraph, nodeMapping, directedEdges, extendNodeCnt, u);
		} else if (parameters["extendStrategy"] == "edgeScores") {
			extend_edgeScores(egoGraph, nodeMapping, directedEdges, extendNodeCnt, u);
		} else {
			throw std::runtime_error("Missing strategy to extend ego-net");
		}

		// Filter edges
		double discardThreshold = std::stod(parameters["discardThreshold"]);
		egoGraph.forEdges([&](node u, node v, edgeweight weight) {
			if (weight < discardThreshold) {
				egoGraph.removeEdge(u, v);
			}
		});

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
		timings["1f)    Copy EgoNet"] += timer.elapsedMicroseconds();


		DEBUG("Cluster EgoNet");
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
		timer.stop();
		timings["1c)    Cluster EgoNet"] += timer.elapsedMicroseconds();

		timer.start();
		egoPartition.compact();
		timer.stop();
		timings["1d)    Compact EgoNet"] += timer.elapsedMicroseconds();


		DEBUG("Build EgoNet");
		timer.start();
		// Insert nodes into ego-net data structure
		for (node i : egoGraph.nodes()) {
			egoNetPartitions[u].emplace(nodeMapping.global(i), egoPartition.subsetOf(i));
		}
		timer.stop();
		timings["1e)    EgoNet subsets"] += timer.elapsedMicroseconds();

		timer.start();
		egoNetPartitionCounts[u] = egoPartition.numberOfSubsets();
		timer.stop();
		timings["1f)    EgoNet subsetCnt"] += timer.elapsedMicroseconds();


		DEBUG("Clean up");
		timer.start();
		nodeMapping.reset();
		timer.stop();
		timings["1g)    Clean up"] += timer.elapsedMicroseconds();
	});
}

void EgoSplitting::extend_edgeScores(Graph &egoGraph,
									 NodeMapping &nodeMapping,
									 const NetworKit::AdjacencyArray &directedEdges,
									 count extendNodeCnt,
									 node u) {
	std::vector<node> secondNeighbors;
	std::vector<std::vector<double>> edgeScores(G.upperNodeIdBound());
	std::vector<node> directNeighbors = nodeMapping.globalNodes();

	// Search for all edges to neighbors of neighbors
	auto isDirectNeighbor = [&](node x) {
		return nodeMapping.isMapped(x);
	};
	for (node v : directNeighbors) {
		edgeScoreGraph.forEdgesOf(v, [&](node, node w, edgeweight weight) {
			if (!isDirectNeighbor(w) && w != u) { // TODO: Adding the node u to the ego-net seems to improve the result
				if (edgeScores[w].empty())
					secondNeighbors.push_back(w);
				edgeScores[w].push_back(weight);
			}
		});
	}

	// Sort nodes by score
	std::string scoreStrategy = parameters["scoreStrategy"];
	auto getScore = [scoreStrategy](std::vector<double> scores) {
		assert(!scores.empty());
		if (scores.size() < 2)
			return 0.0;
		count cnt = scores.size();
		double score = 0.0;
		if (scoreStrategy == "count")
			return static_cast<double>(scores.size());
		if (scoreStrategy == "geometric") {
			std::sort(scores.begin(), scores.end(), std::greater<double>());
			score = (scores[0] + scores[1]) / 4;
			for (count i = 2; i < cnt; ++i) {
				score += scores[i] / std::pow(2, i);
			}
			return score;
		}
		if (scoreStrategy == "numScores") {
			count num_scores = 4;
			for (count i = 0; i < num_scores; ++i) {
				if (cnt > i)
					score += scores[i];
			}
			return score / num_scores;
		}
		throw std::runtime_error("No valid score strategy provided");
	};
	std::vector<std::pair<node, double>> nodeScores;
	for (node v : secondNeighbors) {
		nodeScores.emplace_back(v, getScore(edgeScores[v]));
	}
	std::sort(nodeScores.begin(), nodeScores.end(),
			[](std::pair<node, double> a, std::pair<node, double> b) {
		return a.second > b.second;
	});
	// Remove nodes with score zero
	for (size_t i = 0; i < nodeScores.size(); ++i) {
		if (nodeScores[i].second <= 0) {
			nodeScores.resize(i);
			break;
		}
	}

	// Add nodes to egoGraph
	if (nodeScores.size() > extendNodeCnt)
		nodeScores.resize(extendNodeCnt);
	for (auto pair : nodeScores) {
		nodeMapping.addNode(pair.first);
		egoGraph.addNode();
	}

	// Add edges to egoGraph
	count directNeighborsCnt = directNeighbors.size();
	for (node v : nodeMapping.globalNodes()) {
		edgeScoreGraph.forEdgesOf(v, [&](node, node w, edgeweight weight) {
			if (nodeMapping.isMapped(w)) {
				// Discard edges between neighbors of neighbors
				if (parameters["edgesBetweenNeigNeig"] == "No"
					&& nodeMapping.local(v) >= directNeighborsCnt
					&& nodeMapping.local(w) >= directNeighborsCnt) {
					return;
				}
				if (!egoGraph.hasEdge(nodeMapping.local(v), nodeMapping.local(w)))
					egoGraph.addEdge(nodeMapping.local(v), nodeMapping.local(w), weight);
			}
		});
	}
}

void EgoSplitting::extend_simpleNN(Graph &egoGraph,
								   NodeMapping &nodeMapping,
								   AdjacencyArray const &directedEdges,
								   count extendNodeCnt,
								   node u) {
	count directNeighborsCnt = nodeMapping.nodeCount();
	// Add neighbors of neighbors
	std::vector<node> dirNeighb = nodeMapping.globalNodes();
	for (node v : dirNeighb) {
		auto work = [&](node, node w, edgeweight weight) {
			if (w != u) {
				nodeMapping.addNode(w);
				egoGraph.addNode();
			}
		};
		if (parameters["searchNeigOfNeighInDirectedGraph"] == "Yes") {
			directedEdges.forEdgesOf(v, work);
		} else {
			G.forEdgesOf(v, work);
		}
	}

	DEBUG("Find edges");
	// Add edges to the egoGraph
//	Graph protoEgoGraph(idToNode.size());
	for (node v : nodeMapping.globalNodes()) {
		directedEdges.forEdgesOf(v, [&](node, node w, edgeweight weight) {
			if (nodeMapping.isMapped(w)) {
				// Discard edges between neighbors of neighbors
				if (parameters["discardNeigOfNeigEdgesAtFirst"] == "Yes"
					&& nodeMapping.local(v) >= directNeighborsCnt
					&& nodeMapping.local(w) >= directNeighborsCnt) {
					return;
				}
				if (!egoGraph.hasEdge(nodeMapping.local(v), nodeMapping.local(w)))
					egoGraph.addEdge(nodeMapping.local(v), nodeMapping.local(w), weight);
			}
		});
	}

	Graph extendedGraph = egoGraph;
	if (parameters["discardNonTriangle"] == "Yes") {
		// Discard edges that are not part of a triangle
		extendedGraph = Graph(egoGraph.upperNodeIdBound());
		for (node v : extendedGraph.nodes()) {
			if (!egoGraph.hasNode(v))
				extendedGraph.removeNode(v);
		}

		AdjacencyArray directedEgoGraph(egoGraph);
		std::vector<double> neighborEdgeWeight(egoGraph.upperNodeIdBound(), 0.0);
		egoGraph.forNodes([&](node v) {
			// Mark neighbors
			directedEgoGraph.forEdgesOf(v, [&](node, node w, edgeweight weight) {
				neighborEdgeWeight[w] = weight;
			});

			directedEgoGraph.forEdgesOf(v, [&](node, node w, edgeweight weight) {
				// Always add edges between direct neighbors of u
				if (v < directNeighborsCnt && w < directNeighborsCnt)
					extendedGraph.addEdge(v, w, weight);
				directedEgoGraph.forEdgesOf(w, [&](node, node x, edgeweight weight2) {
					if (neighborEdgeWeight[x] != 0.0) {
						// we have found a triangle v-w-x
						auto tryAddEdge = [&](node u1, node u2, edgeweight edgeWeight) {
							if ((u1 >= directNeighborsCnt || u2 >= directNeighborsCnt) && !extendedGraph.hasEdge(u1, u2))
								extendedGraph.addEdge(u1, u2, edgeWeight);
						};
						tryAddEdge(v, w, weight);
						tryAddEdge(w, x, weight2);
						tryAddEdge(v, x, neighborEdgeWeight[x]);
					}
				});
			});

			// Reset marked neighbors
			directedEgoGraph.forEdgesOf(v, [&](node, node w, edgeweight weight) {
				neighborEdgeWeight[w] = 0.0;
			});
		});
	}

	if (parameters["removeLowDegreeNeig"] == "Yes") {
		// Remove neighbors of neighbors with low degree
		count addedNodes = extendedGraph.nodes().size() - directNeighborsCnt;
		count initMinDegree = count(std::stoi(parameters["minDegreeCleaning"]));
		for (count minDegree = initMinDegree; minDegree <= 10; ++minDegree) {
			auto egoNodes = extendedGraph.nodes();
			for (node v : egoNodes) {
				if (addedNodes <= extendNodeCnt && minDegree > initMinDegree)
					break;
				// Filtering high degree nodes seems to produce worse results, even though these
				// nodes are often part of multiple communities
				if (v >= directNeighborsCnt && extendedGraph.degree(v) < minDegree) {
					extendedGraph.removeNode(v);
					--addedNodes;
				}
			}
		}
	}

	if (parameters["edgesBetweenNeigNeig"] == "Yes") {
		// Add edges between neighbors of neighbors
		for (node v : extendedGraph.nodes()) {
			G.forEdgesOf(nodeMapping.global(v), [&](node, node w, edgeweight weight) {
				node id_w = nodeMapping.local(w);
				if (extendedGraph.hasNode(id_w) && v >= directNeighborsCnt
					&& id_w >= directNeighborsCnt && !extendedGraph.hasEdge(v, id_w)) {
					extendedGraph.addEdge(v, id_w, weight);
				}
			});
		}
	}

	egoGraph = extendedGraph;
}

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

void EgoSplitting::setParameters(std::map<std::string, std::string> parameters) {
	for (auto &x : parameters) {
		this->parameters[x.first] = x.second;
	}
}

} /* namespace NetworKit */
