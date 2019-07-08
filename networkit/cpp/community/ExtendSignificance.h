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
#include "../structures/AdjacencyArray.h"
#include "../structures/Partition.h"
#include "../base/Algorithm.h"
#include "EgoSplitting.h"
#include "ExtendScore.h"
#include "../structures/MemoizationTable.h"

namespace NetworKit {

struct GroupStubs {
	GroupStubs() = default;
	explicit GroupStubs(int size) : groupTotal(size), groupOutgoing(size), externalStubs(size),
	                                externalNodes(size) {
	}

	std::vector<count> groupTotal;
	std::vector<count> groupOutgoing;
	std::vector<count> externalStubs;
	std::vector<count> externalNodes;
};

class ExtendSignificance : public ExtendScore {

public:
	ExtendSignificance(const EgoNetData &egoNetData,
	                   const Partition &basePartition);

	void run() override;

	std::string toString() const override;

	bool isParallel() const override;

private:
	const Partition &basePartition;
	MemoizationTable<double> &sigTable;
	Graph coarseGraph;
	std::vector<std::vector<count>> edgesToGroups; // For each candidate: number of edges to the groups
	std::map<node, std::vector<node> > coarseToEgo;
	std::vector<node> egoToCoarse;
	std::vector<count> coarseSizes;
	int orderStatPos = 0;
	count numGroups = 0;
	std::vector<std::pair<count, node>> candidatesSorted;
	double scorePenalty = 0.0;
	GroupStubs groupStubs;
	std::vector<node> significantGroup;
	std::unordered_set<node> addedCandidates;
	// Algorithm parameters
	const bool useSigMemo;
	const bool mergeGroups;
	const bool sortGroupsStrat;
	const double maxSignificance;
	const count maxGroupCnt;
	const count minEdgesToGroup;

	void
	checkCandidates(const std::string &t_prefix);

	std::vector<node> getCandidates();

	node addExtEdges(std::vector<node> &candidates);

	std::vector<std::pair<count, node>> sortCandidatesByEdges(const std::vector<node> &candidates) const;

	void secondRound();

	void createCoarseGraph();

	double
	calcScore(int nodeDegree, int kIn, int grOut, int groupExtStubs, int extNodes) const;

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

	GroupStubs calcGroupStubsCounts() const;

	std::vector<std::pair<double, node>> sortGroupsByEdges(node v) const;

	void
	checkCandidate(const std::string &t_prefix, Aux::Timer &timer, node v);

	void updateCandidates();
};

} /* namespace NetworKit */

#endif //EXTENDSIGNIFICANCE_H
