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

#include "../structures/UnionFind.h"
#include "../structures/Partition.h"
#include "../structures/AdjacencyArray.h"
#include "../components/ConnectedComponents.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/SignalHandling.h"
#include "../auxiliary/Timer.h"
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
	parameters["discardNeigOfNeigEdgesAtFirst"] = "No";
	parameters["discardNonTriangle"] = "No";
	parameters["minDegreeCleaning"] = "4";
	parameters["maxDegreeCleaning"] = "999999999";
	parameters["addNodesFactor"] = "0.5";
	parameters["removeLowDegreeNeig"] = "No";
	parameters["edgesBetweenNeigNeig"] = "No";
	parameters["extendEgoNet"] = "No";
	parameters["searchNeigOfNeighInDirectedGraph"] = "No";
	parameters["searchEdgesInDirected"] = "No";
}

void EgoSplitting::run() {
	Aux::SignalHandler handler;
	Aux::Timer timer;
	timings.clear();

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

void EgoSplitting::createEgoNets() {
	AdjacencyArray directedEdges(G); // store each undirected edge as one directed edge
	// Assign IDs to the neighbours
	std::vector<count> nodeToId(G.upperNodeIdBound(), none);

	Aux::SignalHandler handler;
	Aux::Timer timer;

	G.forNodes([&](node u) {
		handler.assureRunning();
		DEBUG("Create EgoNet for Node ", u, "/", G.upperNodeIdBound());
		timer.start();
		count degree = G.degree(u);
		// Assign IDs from 0 to degree-1 to nodes of the ego-net
		std::vector<node> idToNode;
		auto addNode = [&](node x){
			if (nodeToId[x] == none) {
				nodeToId[x] = idToNode.size();
				idToNode.push_back(x);
			}
		};
		// Add neighbors
		G.forEdgesOf(u, [&](node, node v) {
			addNode(v);
		});
		if (parameters["extendEgoNet"] == "Yes") {
			// Add neighbors of neighbors
			std::vector<node> dirNeighb = idToNode;
			for (node v : dirNeighb) {
				auto work = [&](node, node w, edgeweight weight) {
					addNode(w);
				};
				if (parameters["searchNeigOfNeighInDirectedGraph"] == "Yes") {
					directedEdges.forEdgesOf(v, work);
				} else {
					G.forEdgesOf(v, work);
				}
			}
		}
		timer.stop();
		timings["1a)    Find nodes"] += timer.elapsedMicroseconds();


		DEBUG("Find edges");
		timer.start();
		// Add edges to the egoGraph
		Graph protoEgoGraph(idToNode.size());
		for (node v : idToNode) {
			auto work = [&](node, node w, edgeweight weight) {
				if (nodeToId[w] != none) {
					// Discard edges between neighbors of neighbors
					if (parameters["discardNeigOfNeigEdgesAtFirst"] == "Yes"
						&& nodeToId[v] >= degree && nodeToId[w] >= degree) {
						return;
					}
					if (!protoEgoGraph.hasEdge(nodeToId[v], nodeToId[w]))
						protoEgoGraph.addEdge(nodeToId[v], nodeToId[w], weight);
				}
			};
			if (parameters["searchEdgesInDirected"] == "Yes") {
				directedEdges.forEdgesOf(v, work);
			} else {
				G.forEdgesOf(v, work);
			}
		}
		timer.stop();
		timings["1b)    Find edges"] += timer.elapsedMicroseconds();

		timer.start();
		Graph egoGraph = protoEgoGraph;
		if (parameters["discardNonTriangle"] == "Yes") {
			// Discard edges that are not part of a triangle
			egoGraph = Graph(protoEgoGraph.upperNodeIdBound());
			for (node v : egoGraph.nodes()) {
				if (!protoEgoGraph.hasNode(v))
					egoGraph.removeNode(v);
			}

			AdjacencyArray directedEgoGraph(protoEgoGraph);
			std::vector<double> neighborEdgeWeight(protoEgoGraph.upperNodeIdBound(), 0.0);
			protoEgoGraph.forNodes([&](node v) {
				// Mark neighbors
				directedEgoGraph.forEdgesOf(v, [&](node, node w, edgeweight weight){
					neighborEdgeWeight[w] = weight;
				});

				directedEgoGraph.forEdgesOf(v, [&](node, node w, edgeweight weight){
					// Always add edges between direct neighbors of u
					if (v < degree && w < degree)
						egoGraph.addEdge(v, w, weight);
					directedEgoGraph.forEdgesOf(w, [&](node, node x, edgeweight weight2) {
						if (neighborEdgeWeight[x] != 0.0) {
							// we have found a triangle v-w-x
							auto tryAddEdge = [&](node u1, node u2, edgeweight edgeWeight) {
								if ((u1 >= degree || u2 >= degree) && !egoGraph.hasEdge(u1, u2))
									egoGraph.addEdge(u1, u2, edgeWeight);
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
		timer.stop();
		timings["1c)    Discard non-triangle edges"] += timer.elapsedMicroseconds();


		timer.start();
		if (parameters["removeLowDegreeNeig"] == "Yes") {
			// Remove neighbors of neighbors with low degree
			count addedNodes = egoGraph.nodes().size() - degree;
			assert(std::stoi(parameters["minDegreeCleaning"]) > 0);
			count initMinDegree = count(std::stoi(parameters["minDegreeCleaning"]));
			assert(std::stoi(parameters["maxDegreeCleaning"]) > 0);
			count maxDegree = count(std::stoi(parameters["maxDegreeCleaning"]));
			assert(std::stod(parameters["addNodesFactor"]) > 0.0);
			double addNodesFactor = std::stod(parameters["addNodesFactor"]);
			for (count minDegree = initMinDegree; minDegree <= 10; ++minDegree) {
				auto egoNodes = egoGraph.nodes();
				for (node v : egoNodes) {
					if (addedNodes <= addNodesFactor * degree && minDegree > initMinDegree)
						break;
					// Filtering high degree nodes seems to produce worse results, even though these
					// nodes are often part of multiple communities
					if (v >= degree &&
						(egoGraph.degree(v) < minDegree || egoGraph.degree(v) >= maxDegree)) {
						egoGraph.removeNode(v);
						--addedNodes;
					}
				}
			}
		}
		timer.stop();
		timings["1d)    Remove low degree neighbors"] += timer.elapsedMicroseconds();

		timer.start();
		if (parameters["edgesBetweenNeigNeig"] == "Yes") {
			// Add edges between neighbors of neighbors
			for (node v : egoGraph.nodes()) {
				G.forEdgesOf(idToNode[v], [&](node, node w, edgeweight weight) {
					node id_w = nodeToId[w];
					if (egoGraph.hasNode(id_w)) {
						if (v >= degree && id_w >= degree &&!egoGraph.hasEdge(v, nodeToId[w]))
							egoGraph.addEdge(v, nodeToId[w], weight);
					}
				});
			}
		}
		timer.stop();
		timings["1e)    Add edges between neigbors of neigbors"] += timer.elapsedMicroseconds();

		timer.start();
		// Get EgoNet with gloabl node ids
		egoNets[u] = Graph(G.upperNodeIdBound());
		egoGraph.forEdges([&](node v, node w){
			egoNets[u].addEdge(idToNode[v], idToNode[w]);
		});
		for (node v : egoNets[u].nodes()) {
			if (!egoGraph.hasNode(nodeToId[v]))
				egoNets[u].removeNode(v);
		}
		timer.stop();
		timings["1f)    Copy EgoNet"] += timer.elapsedMicroseconds();


		DEBUG("Cluster EgoNet");
		timer.start();
		// Cluster ego-net with the local cluster algorithm
		Partition egoPartition;
		if (egoGraph.numberOfEdges() > 0) {
			Graph egoCopy(egoGraph);
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
			egoNetPartitions[u].emplace(idToNode[i], egoPartition.subsetOf(i));
		}
		timer.stop();
		timings["1e)    EgoNet subsets"] += timer.elapsedMicroseconds();

		timer.start();
		egoNetPartitionCounts[u] = egoPartition.numberOfSubsets();
		timer.stop();
		timings["1f)    EgoNet subsetCnt"] += timer.elapsedMicroseconds();


		DEBUG("Clean up");
		timer.start();
		// Reset IDs
		for (node v : idToNode) {
			nodeToId[v] = none;
		}
		for (node i : nodeToId)
			assert(i == none);
		timer.stop();
		timings["1g)    Clean up"] += timer.elapsedMicroseconds();
	});
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
	auto iso2 = numIsolatedNodes(personaGraph);
	assert(iso2 == 0);
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
