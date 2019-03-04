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

EgoSplitting::EgoSplitting(const Graph &G,
						   std::function<Partition(Graph &)> localClusterAlgo,
						   std::function<Partition(Graph &)> globalClusterAlgo,
						   const Cover &groundTruth)
		: G(G),
		  localClusterAlgo(std::move(localClusterAlgo)),
		  globalClusterAlgo(std::move(globalClusterAlgo)),
		  groundTruth(groundTruth) {
	init();
}

void EgoSplitting::init() {
	egoNetPartitions.resize(G.upperNodeIdBound());
	personaOffsets.resize(G.upperNodeIdBound() + 1, 0);
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

	// Store number of partitions of the ego-net
	for (std::size_t i = 0; i < G.upperNodeIdBound(); ++i) {
		egoNetPartitions[i].emplace(none, 0);
	}

	Aux::SignalHandler handler;
	Aux::Timer timer;

	count allPartitions = 0;
	count partitions = 0;
	count allComponents = 0;
	count components = 0;
	std::pair<count, double> f1Score(0, 0.0);
	egoPartitions.resize(G.upperNodeIdBound());

	G.forNodes([&](node u) {
		handler.assureRunning();
		DEBUG("Create EgoNet for Node ", u, "/", G.upperNodeIdBound());
		timer.start();
		count degree = G.degree(u);
		// Assign IDs from 0 to degree-1 to neighbors
		std::vector<node> idToNode(degree);
		{
			index i = 0;
			G.forEdgesOf(u, [&](node, node v) {
				idToNode[i] = v;
				nodeToId[v] = i++;
			});
			assert(i == degree);
		}
		timer.stop();
		timings["1a)    Assign IDs"] += timer.elapsedMicroseconds();


		DEBUG("Find triangles");
		UnionFind unionFind(degree);
		timer.start();
		// Find all triangles and add the edges to the egoGraph
		Graph egoGraph(degree);
		G.forEdgesOf(u, [&](node, node v) {
			directedEdges.forEdgesOf(v, [&](node, node w, edgeweight weight) {
				if (nodeToId[w] != none) {
					// we have found a triangle u-v-w
					egoGraph.addEdge(nodeToId[v], nodeToId[w], weight);
					unionFind.merge(nodeToId[v], nodeToId[w]);
				}
			});
		});
		timer.stop();
		timings["1b)    Find triangles"] += timer.elapsedMicroseconds();


		DEBUG("Cluster EgoNet");
		timer.start();
		// Cluster ego-net with the local cluster algorithm
		Partition egoPartition;
		if (egoGraph.numberOfEdges() > 0) {
			Graph egoCopy(egoGraph);
			egoPartition = localClusterAlgo(egoCopy); // Python moves the graph so we need a copy
		} else {
			egoPartition = unionFind.toPartition();
		}
		timer.stop();
		timings["1c)    Cluster EgoNet"] += timer.elapsedMicroseconds();

		timer.start();
		egoPartition.compact();
		timer.stop();
		timings["1d)    Compact EgoNet"] += timer.elapsedMicroseconds();

		// Count partitions of EgoNet
		for (count size : egoPartition.subsetSizes()) {
			if (size > 1)
				++partitions;
			++allPartitions;
		}
		auto comp = unionFind.toPartition();
		for (count size : comp.subsetSizes()) {
			if (size > 1)
				++components;
			++allComponents;
		}

		// Compare partitions with ground truth
		// For each partition in the EgoNet, check if a community exists in the ground truth
		if (groundTruth.upperBound() > 1 && egoGraph.numberOfNodes() > 1) {
			Cover egoGroundTruth(degree);
			egoGroundTruth.setUpperBound(groundTruth.upperBound());
			egoGraph.forNodes([&](node u) {
				auto globalNode = idToNode[u];
				auto communites = groundTruth.subsetsOf(globalNode);
				for (index c : communites) {
					egoGroundTruth.addToSubset(c, u);
				}
			});
			Cover egoCover = Cover(egoPartition);
		}


		DEBUG("Build EgoNet");
		timer.start();
		// Insert nodes into ego-net data structure
		for (index i = 0; i < degree; ++i) {
			egoNetPartitions[u].emplace(idToNode[i], egoPartition.subsetOf(i));
		}
		timer.stop();
		timings["1e)    EgoNet subsets"] += timer.elapsedMicroseconds();

		timer.start();
		egoNetPartitions[u][none] = egoPartition.numberOfSubsets();
		timer.stop();
		timings["1f)    EgoNet subsetCnt"] += timer.elapsedMicroseconds();


		DEBUG("Clean up");
		timer.start();
		// Reset IDs
		for (node v : idToNode) {
			nodeToId[v] = none;
		}
		timer.stop();
		timings["1g)    Clean up"] += timer.elapsedMicroseconds();
	});
	executionInfo["egoF1Score"] = f1Score.second;
	executionInfo["allPartitions"] = allPartitions / (double) G.numberOfNodes();
	executionInfo["twoPlusPartitions"] = partitions / (double) G.numberOfNodes();
	executionInfo["allComponents"] = allComponents / (double) G.numberOfNodes();
	executionInfo["twoPlusComponents"] = components / (double) G.numberOfNodes();
}

void EgoSplitting::splitIntoPersonas() {
	count sum = 0;
	for (index i = 0; i < G.upperNodeIdBound(); ++i) {
		personaOffsets[i] = sum;
		sum += egoNetPartitions[i][none];
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

} /* namespace NetworKit */
