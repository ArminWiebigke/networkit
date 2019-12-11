/*
 * ExtendSignificance.h
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#ifndef EXTENDSIGNIFICANCE_H
#define EXTENDSIGNIFICANCE_H

#include <vector>

#include "../../graph/Graph.h"
#include "../../structures/NodeMapping.h"
#include "../../auxiliary/Timings.h"
#include "../../structures/LowToHighDirectedGraph.h"
#include "../../structures/Partition.h"
#include "../../base/Algorithm.h"
#include "EgoSplitting.h"
#include "ExtendEgoNetStrategy.h"
#include "../../structures/MemoizationTable.h"
#include "../cleanup/SignificanceCalculator.h"

namespace NetworKit {

struct GroupProperties {
	count groupTotal;
	count groupOutgoing;
	count externalStubs;
	count externalNodes;
};

class ExtendSignificance : public ExtendEgoNetStrategy {

public:
	ExtendSignificance(EgoNetData &egoNetData,
	                   const Partition &basePartition, count maxCandidates,
	                   const Graph &egoGraph, node egoNode);

	void run() override;

	std::string toString() const override;

	bool isParallel() const override;

private:
	const Partition &basePartition;
	MemoizationTable<double> &sigTable;
	SignificanceCalculator &significanceCalculator;
	Graph coarseGraph;
	SparseVector<std::vector<count>> &edgesToGroups; // For each candidate: number of edges to the groups
	std::map<node, std::vector<node> > coarseToEgo;
	std::vector<node> egoToCoarse;
	std::vector<count> coarseSizes;
	count numGroups = 0;
	std::vector<node> candidatesSorted;
	std::vector<GroupProperties> groupProperties;
	SparseVector<node> &significantGroup;
	std::unordered_set<node> addedCandidates;
	// Algorithm parameters
	const bool useSigMemo;
	const bool mergeGroups;
	const bool sortGroupsBySignificance;
	const double maxSignificance;
	const count maxGroupCnt;
	const count minEdgesToGroup;  // Minimum number of edges to check for significance
	const bool onlyCheckMaxCandidates;

	void
	findSignificantCandidates(const std::string &t_prefix);

	std::vector<node> getCandidates();

	void insertOutgoingEdgesIntoCoarseGraph(std::vector<node> &candidates);

	std::vector<node> sortCandidatesByEdges(const std::vector<node> &candidates) const;

	void createCoarseGraph();

	double
	calculateSScore(count nodeDegree, count kIn, count grOut, count groupExtStubs, count extNodes) const;

	bool
	checkSignificanceToMergedGroups(const std::string &t_prefix, Aux::Timer &timer, node candidate,
	                                std::vector<std::pair<double, node>> &numEdgesToGroups,
	                                const std::vector<double> &groupSigs);

	bool addIfSignificant(node v, double significance, node group);

	bool
	checkSignificanceToSingleGroups(node candidate,
	                                const std::vector<std::pair<double, node>> &numEdgesToGroups,
	                                std::vector<double> &significanceToGroups);

	std::vector<GroupProperties> calculateGroupProperties() const;

	std::vector<std::pair<double, node>> sortGroupsByEdges(node candidate) const;

	void
	checkCandidate(const std::string &t_prefix, node candidate);

	void updateCandidates();

	bool enoughSignificantCandidates() const;

	void resetData();

	void setMemoizationFunction();
};

} /* namespace NetworKit */

#endif //EXTENDSIGNIFICANCE_H
