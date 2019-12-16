/*
 * ExtendByScore.cpp
 *
 * Created: 2019-06-19
 * Author: Armin Wiebigke
 */

#include <networkit/community/egosplitting/ExtendByScore.hpp>

namespace NetworKit {

ExtendByScore::ExtendByScore(EgoNetData &egoNetData, count maxCandidates,
                             const Graph &egoGraph, node egoNode)
        : ExtendEgoNetStrategy(egoNetData, maxCandidates, egoGraph, egoNode),
          nodeScores(egoNetData.nodeScores), scoreStrategy(parameters.at("Edges Score Strategy")),
          significanceCalculator(egoNetData.significanceCalculator) {
    if (nodeScores.upperBound() < G.upperNodeIdBound())
        nodeScores.setUpperBound(G.upperNodeIdBound());
}

void ExtendByScore::run() {
    searchForCandidates();

    std::vector<NodeAndScore> candidatesAndScores = calculateScores();

    takeBestCandidates(candidatesAndScores);

    nodeScores.reset();

    hasRun = true;
}

void ExtendByScore::searchForCandidates() {
    const std::vector<node>& egoNetNodes = egoMapping.globalNodes();
    auto isInEgoNet = [&](node x) {
        return egoMapping.isMapped(x);
    };
    count internalStubs = 0;
    outgoingStubs = 0;
    auto countEdges = [&](node, node neighbor, edgeweight weight) {
        if (isInEgoNet(neighbor)) {
            ++internalStubs;
        } else {
            ++outgoingStubs;

            // Exclude ego-net nodes and ego-node as candidates
            if (neighbor != egoNode) {
                if (!nodeScores.indexIsUsed(neighbor)) {
                    nodeScores.insert(neighbor, 0.0);
                }
                nodeScores[neighbor] += weight;
            }
        }
    };
    for (node egoNetNode : egoNetNodes) {
        G.forEdgesOf(egoNetNode, countEdges);
    }
    externalStubs = G.numberOfEdges() * 2 - internalStubs - outgoingStubs;
}

template<typename F>
std::vector<ExtendByScore::NodeAndScore>
ExtendByScore::calculateScoresImpl(F calculateScore) const {
    std::vector<NodeAndScore> candidatesAndScores;
    candidatesAndScores.reserve(nodeScores.size());
    for (node candidate : nodeScores.insertedIndexes()) {
        double numEdges = nodeScores[candidate];
        if (numEdges >= 3) {
            double score = calculateScore(candidate, numEdges);
            candidatesAndScores.emplace_back(candidate, score);
        }
    }
    return candidatesAndScores;
}

std::vector<ExtendByScore::NodeAndScore>
ExtendByScore::calculateScores() const {
    std::vector<NodeAndScore> candidatesAndScores;
    if (scoreStrategy == "Edges pow 2 div Degree") {
        candidatesAndScores = calculateScoresImpl([&](node v, count numEdges) {
            return numEdges * numEdges * 1.0 / G.degree(v);
        });
    } else if (scoreStrategy == "constant") {
        candidatesAndScores = calculateScoresImpl([](node, count) { return 1.0; });
    } else if (scoreStrategy == "Random") {
        candidatesAndScores = calculateScoresImpl([](node, count) { return Aux::Random::real(); });
    } else if (scoreStrategy == "Edges" || scoreStrategy == "none") {
        candidatesAndScores = calculateScoresImpl([](node, count numEdges) { return numEdges; });
    } else if (scoreStrategy == "Edges div Degree") {
        candidatesAndScores = calculateScoresImpl([&](node v, count numEdges) {
            return numEdges * 1.0 / G.degree(v);
        });
    } else if (scoreStrategy == "Significance") {
        candidatesAndScores = calculateScoresImpl([&](node v, count numEdges) {
            double rScore = significanceCalculator.rScore(G.degree(v), numEdges, outgoingStubs, externalStubs);
            return -rScore; // low r-score is better
        });
    } else {
        throw std::runtime_error(scoreStrategy + " is not a valid score strategy!");
    }
    return candidatesAndScores;
}

void ExtendByScore::takeBestCandidates(std::vector<NodeAndScore> &candidatesAndScores) {
    if (candidatesAndScores.size() > maxExtendedNodes) {
        std::nth_element(candidatesAndScores.begin(),
                 candidatesAndScores.begin() + maxExtendedNodes,
                 candidatesAndScores.end(),
            [](NodeAndScore a, NodeAndScore b) {
                return a.second > b.second;
            });
        candidatesAndScores.resize(maxExtendedNodes);
    }
    for (NodeAndScore nodeAndScore : candidatesAndScores) {
        significantCandidates.push_back(nodeAndScore.first);
    }
}

bool ExtendByScore::isParallel() const {
    return false;
}

std::string ExtendByScore::toString() const {
    return "ExtendByScore";
}

} /* namespace NetworKit */
