/*
 * EgoSplitting.cpp
 *
 * Created: 2019-05-28
 * Author: Armin Wiebigke
 */

#include <unordered_map>

#include "SLPA.h"
#include "../auxiliary/Random.h"

namespace NetworKit {


SLPA::SLPA(const Graph &graph, double threshold, count numIterations)
		: Algorithm(), graph(graph), numIterations(numIterations), threshold(threshold),
		  cover(graph.upperNodeIdBound()), partition(graph.upperNodeIdBound()),
		  nodesMemory(graph.upperNodeIdBound()) {
	for (node u : graph.nodes())
		addLabelTo(u, u);
}

SLPA::SLPA(const Graph &graph, const Partition &basePartition, double threshold,
           count numIterations)
		: Algorithm(), graph(graph), numIterations(numIterations), threshold(threshold),
		  cover(graph.upperNodeIdBound()), partition(graph.upperNodeIdBound()),
		  nodesMemory(graph.upperNodeIdBound()) {
	for (node u : graph.nodes())
		addLabelTo(u, basePartition.subsetOf(u));
}

void SLPA::run() {
	for (count it = 0; it < numIterations; ++it) {
		graph.forNodesInRandomOrder([&](node u) {
			auto labelCnts = listen(u);
			if (labelCnts.empty())
				return;
			label selectedLabel = selectLabel(labelCnts);
			addLabelTo(u, selectedLabel);
		});
	}
	createResult();
	hasRun = true;
}

std::string SLPA::toString() const {
	return "SLPA";
}


bool SLPA::isParallel() const {
	return false;
}

Cover SLPA::getCover() {
	return cover;
}

Partition SLPA::getPartition() {
	return partition;
}

LabelCounts SLPA::listen(node u) {
	LabelCounts labelCnts;
	graph.forNeighborsOf(u, [&](node neighbor) {
		label randLabel = sendLabel(neighbor);
		labelCnts[randLabel] += 1;
	});
	return labelCnts;
}

// TODO: Select randomly
label SLPA::selectLabel(const LabelCounts& labelCounts) {
//	auto greaterCount = [](std::pair<label, count> a, std::pair<label, count> b) {
//		return a.second < b.second;
//	};
//	auto maxLabel = std::max_element(labelCounts.begin(), labelCounts.end(),
//	                                 greaterCount);
//	return maxLabel->first;
	std::vector<label> labels;
	std::vector<count> weights;
	for (auto p : labelCounts) {
		labels.push_back(p.first);
		weights.push_back(p.second);
	}
	std::discrete_distribution<index> dist(weights.begin(), weights.end());
	index randIndex = dist(Aux::Random::getURNG());
	label randomLabel = labels[randIndex];
	return randomLabel;
}

void SLPA::createResult() {
	cover.setUpperBound(graph.upperNodeIdBound());
	partition.setUpperBound(graph.upperNodeIdBound());
	graph.forNodes([&](node u) {
		LabelCounts labelCnts = nodesMemory[u];
		count labelCntSum = 0;
		for (auto it : labelCnts)
			labelCntSum += it.second;
		label bestLabel = none;
		count bestCnt = 0;
		for (auto it : labelCnts) {
			label l = it.first; // TODO: in one assignment
			count cnt = it.second;
			if (static_cast<double>(cnt) / labelCntSum >= threshold) {
				cover.addToSubset(l, u);
				if (cnt > bestCnt) {
					bestCnt = cnt;
					bestLabel = l;
				}
			}
		}
		if (bestLabel != none)
			partition.addToSubset(bestLabel, u);
		else
			partition.toSingleton(u);
	});
}

label SLPA::sendLabel(node u) const {
	std::vector<label> labels;
	std::vector<count> weights;
	for (auto p : nodesMemory[u]) {
		labels.push_back(p.first);
		weights.push_back(p.second);
	}
	std::discrete_distribution<index> dist(weights.begin(), weights.end());
	index randIndex = dist(Aux::Random::getURNG());
	label randomLabel = labels[randIndex];
	return randomLabel;
}

void SLPA::addLabelTo(node u, label l) {
	nodesMemory[u][l] += 1;
}


} /* namespace NetworKit */

