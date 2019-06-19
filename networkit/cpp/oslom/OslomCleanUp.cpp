/*
 * OslomCleanUp.cpp
 *
 * Created: 2019-06-06
 * Author: Armin Wiebigke
 */

#include "OslomCleanUp.h"
#include "Parameters.h"
#include "Stochastics.h"
#include "LogTable.h"
#include "OslomNetGlobal.h"


namespace NetworKit {


OslomCleanUp::OslomCleanUp(const Graph &graph, const Cover &cover)
		: OslomCleanUp(graph, cover, {}) {

}

OslomCleanUp::OslomCleanUp(const Graph &graph, const Cover &cover,
                           std::vector<std::string> args)
		: Algorithm(), graph(graph), cover(cover), args(std::move(args)) {

}

void setUpParameters(const std::vector<std::string> &args) {
	Parameters *paras = Parameters::get_instance();
	paras->set(args);
}

void setUpStochastics(int log_max) {
	Stochastics::init(log_max);
}

std::map<int, std::map<int, std::pair<int, double> > >
read_networkit_graph(const NetworKit::Graph &graph) {
	std::map<int, std::map<int, std::pair<int, double> > > graph_map;
	graph.forNodes([&](NetworKit::node u) {
		std::map<int, std::pair<int, double> > neighbor_map;
		graph.forEdgesOf(u, [&](NetworKit::node, NetworKit::node v, NetworKit::edgeweight weight) {
			neighbor_map.insert(std::make_pair(v, std::make_pair(1, weight)));
		});
		graph_map[u] = neighbor_map;
	});
	return graph_map;
}

std::vector<std::deque<int>> read_networkit_cover(const NetworKit::Cover &cover) {
	std::vector<std::deque<int>> communities(cover.upperBound());
	cover.forEntries([&](NetworKit::node u, std::set<NetworKit::index> cs) {
		for (auto c : cs) {
			communities[c].push_back(u);
		}
	});
	for (size_t i = 0; i < communities.size(); ++i) {
		if (communities[i].empty()) {
			communities[i] = communities.back();
			communities.pop_back();
			--i;
		}
	}
	return communities;
}

void removeDataStructures() {
	Parameters::delete_instance();
	LogFactTable::delete_instance();
}

void OslomCleanUp::run() {
	// Set parameters
	setUpParameters(args);

	// Import graph
	auto graph_map = read_networkit_graph(graph);
	OslomNetGlobal oslom(graph_map);
	setUpStochastics((int) oslom.stubs());

	// Read input cover
	auto communities = read_networkit_cover(cover);

	// Clean up cover
	auto goodModules = oslom.clean_up(communities, graph.upperNodeIdBound());
	oslom.print_statistics(std::cout, goodModules);

	// Create result cover
	resultCover = Cover(graph.upperNodeIdBound());
	resultCover.setUpperBound(cover.upperBound());
	graph.forNodes([&](NetworKit::node u) {
		for (auto c : goodModules.memberships[u])
			resultCover.addToSubset(c, u);
	});

	removeDataStructures();

	hasRun = true;
}

Cover OslomCleanUp::getCover() {
	return resultCover;
}

std::string OslomCleanUp::toString() const {
	return "OslomCleanUp";
}

bool OslomCleanUp::isParallel() const {
	return false;
}


} /* namespace NetworKit */
