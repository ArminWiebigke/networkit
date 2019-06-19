/*
 * ExtendSignificance.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#include <set>

#include "ExtendSignificance.h"
#include "../auxiliary/Timer.h"
#include "../coarsening/ParallelPartitionCoarsening.h"
#include "../oslom/Stochastics.h"
#include "../structures/AdjacencyArray.h"
#include "EgoSplitting.h"

namespace NetworKit {

ExtendSignificance::ExtendSignificance(const EgoNetData &egoNetData,
                                       const Partition &basePartition)
		: ExtendScore(egoNetData), basePartition(basePartition) {

}

void ExtendSignificance::run() {
	Aux::Timer timer;
	timer.start();
	Stochastics::init(G.numberOfEdges());
	ParallelPartitionCoarsening coarsening(egoGraph, basePartition);
	coarsening.run();
	Graph coarseGraph = coarsening.getCoarseGraph();
	auto coarseToEgo = coarsening.getCoarseToFineNodeMapping();
	auto egoToCoarse = coarsening.getFineToCoarseNodeMapping();

	std::vector<count> coarseSizes(coarseGraph.upperNodeIdBound());
	coarseGraph.forNodes([&](node p) {
		coarseSizes[p] = coarseToEgo[p].size();
	});
	addTime(timer, "11s0    Coarsening");

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
		if (!isDirectNeighbor(w) && w != egoNode) {
			assert(!coarseMapping.isMapped(w));
			if (edgeScores[w].empty())
				candidates.insert(w);
			edgeScores[w][egoToCoarse[egoMapping.local(v)]] += weight;
		}
	};

	if (parameters.at("extendOverDirected") == "Yes") {
		for (node v : directNeighbors)
			directedG.forEdgesOf(v, countEdges);
		for (node w : candidates) {
			directedG.forEdgesOf(w, [&](node, node v, edgeweight weight) {
				if (isDirectNeighbor(v))
					edgeScores[w][egoToCoarse[egoMapping.local(v)]] += weight;
			});
		}
	} else {
		for (node v : directNeighbors) {
			G.forEdgesOf(v, countEdges);
		}
	}

//	for (node v : directNeighbors)
//			directedG.forEdgesOf(v, countEdges);
//	for (node w : candidates) {
//		directedG.forEdgesOf(w, [&](node, node v, edgeweight weight) {
//			if (isDirectNeighbor(v))
//				edgeScores[w][egoToCoarse[egoMapping.local(v)]] += weight;
//		});
//	}
//	auto candidatesDir = candidates;
//	candidates.clear();
//	auto edgeScoresDir = edgeScores;
//	edgeScores = std::vector<std::map<node, double>>(G.upperNodeIdBound());
//	for (node v : directNeighbors) {
//		G.forEdgesOf(v, countEdges);
//	}
//	assert(candidatesDir.size() <= candidates.size());
//	for (node w : candidates){
//		for (auto pair : edgeScores[w])
//			assert(edgeScoresDir[w][pair.first] <= edgeScores[w][pair.first]);
//	}

	// Insert candidates into graph
	for (node v : candidates) {
		assert(!coarseMapping.isMapped(v));
		coarseMapping.addNode(v);
		coarseGraph.addNode();
		for (auto p : edgeScores[v]) {
			assert(p.second <= coarseSizes[p.first]);
			coarseGraph.addEdge(coarseMapping.local(v), p.first, p.second);
		}
	}
	addTime(timer, "11s2    Add candidates");

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
	count missingEdges = numEgoEdges - (int) totalWeight;
	coarseGraph.addEdge(externalNode, externalNode, missingEdges);

	// Calculate significances
	result = calcSignficance(externalNode, coarseGraph, coarseMapping, coarseSizes, 0);
	addTime(timer, "11s3a    Calc significane one");

	// Calculate significances, allowing slightly worse candidates (based on the number of already
	// found significant candidates)
	for (auto pair : result) {
		node v = coarseMapping.local(pair.first);
		coarseGraph.increaseWeight(externalNode, externalNode, coarseGraph.weightedDegree(v));
		coarseGraph.removeNode(v);
	}
	double orderedStatPosition = std::stod(parameters.at("orderedStatPos")) * result.size();
	auto scores2 = calcSignficance(externalNode, coarseGraph, coarseMapping, coarseSizes,
	                               orderedStatPosition);
	for (auto pair : scores2)
		result.emplace_back(pair.first,
		                    pair.second - 0.5); // TODO: make sure that this score is low and > 0
	addTime(timer, "11s3b    Calc significane two");

	hasRun = true;
}

std::vector<std::pair<node, double>>
ExtendSignificance::calcSignficance(node externalNode, const Graph &coarseGraph,
                                    const NodeMapping &coarseMapping,
                                    const std::vector<count> &coarseSizes,
                                    double orderedStatPosition) const {
	Aux::Timer timer;
	timer.start();
	// Sort candidates by number of edges
	std::vector<std::pair<count, node>> candidatesSorted;
	candidatesSorted.reserve(coarseGraph.upperNodeIdBound() - externalNode);
	for (node v = externalNode + 1; v < coarseGraph.upperNodeIdBound(); ++v) {
		if (!coarseGraph.hasNode(v))
			continue;
		count edgesIntoEgo = (int) coarseGraph.weightedDegree(v);
		candidatesSorted.emplace_back(edgesIntoEgo, coarseMapping.global(v));
	}
	std::sort(candidatesSorted.rbegin(), candidatesSorted.rend());
	addTime(timer, "11s6    Sort candidates");

	// Discard nodes with less than 3 edges
	for (count i = 0; i < candidatesSorted.size(); ++i) {
		if (candidatesSorted[i].first < 3) {
			candidatesSorted.resize(i);
			break;
		}
	}
	addTime(timer, "11s7    Discard candidates");

	// Calculate stub counts for each group
	std::vector<count> groupTotal(externalNode);
	std::vector<count> groupOutgoing(externalNode);
	std::vector<count> externalStubs(externalNode);
	std::vector<count> externalNodes(externalNode);
	for (count p = 0; p < externalNode; ++p) {
		if (!coarseGraph.hasNode(p))
			continue;
		groupTotal[p] = (int) coarseGraph.weightedDegree(p);
		groupOutgoing[p] = (int) coarseGraph.weightedDegree(p) - 2 * (int) coarseGraph.weight(p, p);
		externalStubs[p] = G.numberOfEdges() * 2 - groupTotal[p];
		externalNodes[p] = G.numberOfNodes() - coarseSizes[p];
	}
	addTime(timer, "11s8    Calculate stub counts");

	index statPosInt = (int) std::floor(orderedStatPosition);
	double statPosWeight = orderedStatPosition - statPosInt;
	if (parameters.at("useSignInterpol") == "No")
		statPosWeight = 0.0;
	double maxSignificance = std::stod(parameters.at("maxSignificance"));
	std::vector<std::pair<node, double>> nodeScores;
	auto calc_score = [&](int nodeDegree, int kIn, int grOut, int extStubs, int extNodes,
	                      int position) {
		// TODO: bootInterval anschauen
		// TODO: Reuse ordered statistics if externalNodes very similar (less than 1% difference)
		assert(grOut <= extStubs);
		if (parameters.at("subtractNodeDegree") == "Yes")
			extStubs -= nodeDegree; // Remove node stubs from the free external stubs
		else
			extStubs -= 1;

		double rScore = Stochastics::compute_simple_fitness(kIn, grOut, extStubs, nodeDegree);
		addTime(timer, "11sc    Calc r-score");
		double orderedStat = Stochastics::order_statistics_left_cumulative(
				extNodes, extNodes - position, rScore);
		addTime(timer, "11sd    Calc ordered statistics");
		return orderedStat;
	};
	auto calcSigScore = [&](int nodeDegree, count kIn, count gOut, count extStubs, count extNodes) {
		double significance;
		if (statPosWeight != 0.0) {
			double s1 = calc_score(nodeDegree, kIn, gOut, extStubs, extNodes, statPosInt);
			double s2 = calc_score(nodeDegree, kIn, gOut, extStubs, extNodes, statPosInt + 1);
			if (s1 == 0.0 || s2 == 0.0)
				significance = 0.0;
			else
				significance = std::exp((1 - statPosWeight) * std::log(s1)
				                        + statPosWeight * std::log(s2));
		} else {
			significance = calc_score(nodeDegree, kIn, gOut, extStubs, extNodes, statPosInt);
		}
		return significance;
	};
	auto addIfSignificant = [&](node v, double significance) {
		if (significance <= maxSignificance) {
			double score = 1 - significance;
			nodeScores.emplace_back(v, score);
			return true;
		}
		return false;
	};

	// Check significance for each candidate
	for (auto pair : candidatesSorted) {
		node v = pair.second;
		addTime(timer, "11sg    Loop");
		count nodeDegree = G.degree(v);
		node vLoc = coarseMapping.local(v);
		addTime(timer, "11s9    Get degree and local");
		// Sort groups by number of edges
		std::vector<std::pair<double, node>> groupEdges;
		auto groups = coarseGraph.neighbors(vLoc);
		for (node p : groups) {
			if (p >= externalNode)
				continue;
			edgeweight numEdges = coarseGraph.weight(vLoc, p);
			groupEdges.emplace_back(numEdges, p);
			assert(numEdges <= coarseSizes[p]);
		}
		std::sort(groupEdges.rbegin(), groupEdges.rend());

		// Check significance to single groups
		std::vector<double> groupSigs(externalNode);
		bool added = false;
		count calcedGroups = 0;
		count maxGroupCnt = std::stoi(parameters.at("maxGroupsConsider"));
		for (auto it : groupEdges) {
			if (calcedGroups++ >= maxGroupCnt)
				break;
			count numEdges = (int) it.first;
			node p = it.second;
			addTime(timer, "11sg    Loop");
			if (numEdges < 3)
				break;

			double significance = calcSigScore(nodeDegree, (int) numEdges, groupOutgoing[p],
			                                   externalStubs[p], externalNodes[p]);
			groupSigs[p] = significance;
			added = addIfSignificant(v, significance);
			addTime(timer, "11sh    Add if significant");
			if (added)
				break;
		}
		if (added)
			continue;

		if (parameters.at("signMerge") == "No")
			continue;

		// Chech significance to merged groups
		if (parameters.at("sortGroups") == "significance") {
			std::sort(groupEdges.begin(), groupEdges.end(), [&](std::pair<double, node> a,
			                                                    std::pair<double, node> b) {
				return (groupSigs[a.second]) > groupSigs[b.second];
			});
		}

		auto it = groupEdges.begin();
		node bestGroup = it->second;
		std::set<node> mergedGroups{bestGroup};
		count totalStubsMerged = groupTotal[bestGroup];
		count outStubsMerged = groupOutgoing[bestGroup];
		count extStubsMerged = externalStubs[bestGroup];
		count extNodesMerged = externalNodes[bestGroup];
		count numEdgesMerged = (int) it->first;

		auto check = [&]() {
			// Check
			count extNodesCheck = G.numberOfNodes();
			for (node w : mergedGroups)
				extNodesCheck -= coarseSizes[w];
			assert(extNodesMerged == extNodesCheck);
			count outStubsCheck = 0;
			count totalStubsCheck = 0;
			for (node w : mergedGroups) {
				coarseGraph.forEdgesOf(w, [&](node, node ext, edgeweight weight) {
					totalStubsCheck += (int) weight;
					if (mergedGroups.count(ext) == 0)
						outStubsCheck += (int) weight;
					if (w == ext)
						totalStubsCheck += (int) weight;
				});
			}
			assert(totalStubsCheck == totalStubsMerged);
			assert(outStubsMerged == outStubsCheck);
			count extStubsCheck = 2 * G.numberOfEdges() - totalStubsCheck;
			assert(extStubsCheck == extStubsMerged);

			assert(extStubsCheck >= outStubsCheck);


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
			assert(g_out == outStubsMerged);
			assert(ext == extStubsMerged);
			return true;
		};
		assert(check());

		for (++it; it < groupEdges.begin() + calcedGroups; ++it) {
			node group = it->second;
			count newInternalStubs = 0;
			for (node w : mergedGroups) {
				newInternalStubs += 2 * (int) coarseGraph.weight(w, group);
			}
			// TODO: Groups have to be connected (?)
			if (newInternalStubs == 0)
				continue;
			// Merge group
			extStubsMerged -= groupTotal[group];
			totalStubsMerged += groupTotal[group];
			outStubsMerged += groupOutgoing[group] - newInternalStubs;
			extNodesMerged -= coarseSizes[group];
			numEdgesMerged += (int) it->first;
			mergedGroups.insert(group);
			addTime(timer, "11sg    Merge groups");
			assert(check());

			// Calculate new significance
			double significance = calcSigScore(nodeDegree, numEdgesMerged, outStubsMerged,
			                                   extStubsMerged, extNodesMerged);
			added = addIfSignificant(v, significance);
			addTime(timer, "11sh    Add if significant");
			if (added) {
//				std::cout << mergedGroups.size() << " merged groups " << nodeScores.back().second
//				          << " (" << numEdgesMerged << "/" << nodeDegree << ") -> "
//				          << G.numberOfNodes() - extNodesMerged
//				          << std::endl;
				break;
			}
		}
		addTime(timer, "11sg    Loop");
	}
	return nodeScores;
}

std::string ExtendSignificance::toString() const {
	return "ExtendSignificance";
}

bool ExtendSignificance::isParallel() const {
	return false;
}

} /* namespace NetworKit */