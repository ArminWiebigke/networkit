/*
 * ExtendSignificance.h
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#ifndef EXTENDSIGNIFICANCE_H
#define EXTENDSIGNIFICANCE_H

#include <vector>

#include "../graph/Graph.h"
#include "../structures/NodeMapping.h"
#include "../auxiliary/Timings.h"
#include "../structures/LowToHighDirectedGraph.h"
#include "../structures/Partition.h"
#include "../base/Algorithm.h"
#include "EgoSplitting.h"
#include "ExtendEgoNetStrategy.h"
#include "../structures/MemoizationTable.h"
#include "../cleanup/StochasticSignificance.h"

namespace NetworKit {

struct GroupProperties {
	GroupProperties() = default;
	explicit GroupProperties(int size) : groupTotal(size), groupOutgoing(size), externalStubs(size),
	                                     externalNodes(size) {
	}

	std::vector<count> groupTotal;
	std::vector<count> groupOutgoing;
	std::vector<count> externalStubs;
	std::vector<count> externalNodes;
};

class ExtendSignificance : public ExtendEgoNetStrategy {

public:
	ExtendSignificance(const EgoNetData &egoNetData,
	                   const Partition &basePartition, count maxCandidates,
	                   const Graph &egoGraph, node egoNode);

	void run() override;

	std::string toString() const override;

	bool isParallel() const override;

private:
	const Partition &basePartition;
	MemoizationTable<double> &sigTable;
	const StochasticSignificance &stochasticSignificance;
	Graph coarseGraph;
	SparseVector<std::vector<count>> &edgesToGroups; // For each candidate: number of edges to the groups
	std::map<node, std::vector<node> > coarseToEgo;
	std::vector<node> egoToCoarse;
	std::vector<count> coarseSizes;
	count numGroups = 0;
	std::vector<node> candidatesSorted;
	GroupProperties groupProperties;
	SparseVector<node> &significantGroup;
	std::unordered_set<node> addedCandidates;
	// Algorithm parameters
	const bool useSigMemo;
	const bool mergeGroups;
	const bool sortGroupsStrat;
	const double maxSignificance;
	const count maxGroupCnt;
	const count minEdgesToGroup;
	const bool onlyCheckMaxCandidates;

	void
	findSignificantCandidates(const std::string &t_prefix);

	std::vector<node> getCandidates();

	node addExtEdges(std::vector<node> &candidates);

	std::vector<node> sortCandidatesByEdges(const std::vector<node> &candidates) const;

	void secondRound();

	void createCoarseGraph();

	double
	calcScore(count nodeDegree, count kIn, count grOut, count groupExtStubs, count extNodes) const;

	bool
	checkMergedGroups(const std::string &t_prefix, Aux::Timer &timer, node v,
	                  std::vector<std::pair<double, node>> &groupEdges,
	                  const std::vector<double> &groupSigs);

	bool addIfSignificant(node v, double significance, node group);

	bool
	checkSingleGroups(node v,
	                  const std::vector<std::pair<double, node>> &groupEdges,
	                  std::vector<double> &groupSigs);

	void removeEgoNodeCandidate(std::vector<node> &candidates);

	GroupProperties calcGroupStubsCounts() const;

	std::vector<std::pair<double, node>> sortGroupsByEdges(node v) const;

	void
	checkCandidate(const std::string &t_prefix, Aux::Timer &timer, node v);

	void updateCandidates();

	bool enoughSignificantCandidates() const;

	void resetData();

	void setMemoizationFunction() const;
};

} /* namespace NetworKit */

#endif //EXTENDSIGNIFICANCE_H
