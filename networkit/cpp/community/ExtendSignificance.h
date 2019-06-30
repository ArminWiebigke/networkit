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
	std::vector<std::vector<double>> edgeScores;
	std::map<node, std::vector<node> > coarseToEgo;
	std::vector<node> egoToCoarse;
	std::vector<count> coarseSizes;
	const bool subtractNodeDegree;
	const bool useSigMemo;
	const bool mergeGroups;
	const bool sortGroupsStrat;
	const double maxSignificance;
	const count maxGroupCnt;
	const count minEdgesToGroup;

	std::vector<std::pair<node, double>>
	calcSignficance(node numGroups, double orderedStatPosition,
	                const std::string &t_prefix,
	                const std::vector<std::pair<count, node>> &candidatesSorted) const;

	std::vector<node> getCandidates();

	node addExtEdges(std::vector<node> &candidates);

	std::vector<std::pair<count, node>> sortCandidates(const std::vector<node> &candidates) const;

	void secondSigRound(node externalNode,
	                    std::vector<std::pair<count, node>> &candidatesSorted,
	                    double orderedStatPosition);

	void createCoarseGraph();

	double
	calcScore(int nodeDegree, int kIn, int grOut, int extStubs, int extNodes, int position) const;

	bool
	checkMergedGroups(const std::string &t_prefix, Aux::Timer &timer,
	                  const GroupStubs &groupStubs, index statPosInt, node v,
	                  std::vector<std::pair<double, node>> &groupEdges,
	                  const std::vector<double> &groupSigs, count calcedGroups,
	                  std::vector<std::pair<node, double>> &nodeScores) const;

	bool addIfSignificant(std::vector<std::pair<node, double>> &nodeScores, node v,
	                      double significance) const;

	bool
	checkSingleGroups(const GroupStubs &groupStubs, index statPosInt,
	                  std::vector<std::pair<node, double>> &nodeScores, node v,
	                  const std::vector<std::pair<double, node>> &groupEdges,
	                  std::vector<double> &groupSigs, count &calcedGroups) const;

	void removeEgoNodeCandidate(std::vector<node> &candidates) const;

	GroupStubs calcGroupStubsCounts(node numGroups) const;

	std::vector<std::pair<double, node>> sortGroupsByEdges(node v) const;

	void
	checkForSignificance(node numGroups, const std::string &t_prefix, Aux::Timer &timer,
	                     const GroupStubs &groupStubs, index statPosInt,
	                     std::vector<std::pair<node, double>> &nodeScores, node v) const;
};

} /* namespace NetworKit */

#endif //EXTENDSIGNIFICANCE_H
