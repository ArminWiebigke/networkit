#include <iostream>
#include <deque>

#include <NetworKit/graph/Graph.h>
#include <NetworKit/Globals.h>
#include <NetworKit/structures/Cover.h>
#include <NetworKit/io/CoverReader.h>
#include <NetworKit/io/EdgeListReader.h>
#include <NetworKit/community/EgoSplitting.h>
#include "OslomNetGlobal.h"
#include "ModuleCollection.h"
#include "LogTable.h"
#include "Stochastics.h"
#include "Parameters.h"

auto read_networkit_graph(const NetworKit::Graph &graph) {
    std::map<int, std::map<int, std::pair<int, double> > > graph_map;
    graph.forNodes([&](NetworKit::node u){
        std::map<int, std::pair<int, double> > neighbor_map;
        graph.forEdgesOf(u, [&](NetworKit::node, NetworKit::node v, NetworKit::edgeweight weight){
            neighbor_map.insert(std::make_pair(v, std::make_pair(1, weight)));
        });
        graph_map[u] = neighbor_map;
    });
    return graph_map;
}

auto read_networkit_cover(const NetworKit::Cover &cover) {
    std::vector<std::deque<int>> communities(cover.upperBound());
    cover.forEntries([&](NetworKit::node u, std::set<NetworKit::index> cs){
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

void set_up_stochastics(int log_max) {
    Stochastics::init();
    LogFactTable *table = LogFactTable::get_instance();
    table->set(log_max);
}

void set_up_parameters(const std::vector<std::string> &args) {
    Parameters *paras = Parameters::get_instance();
    paras->set(args);
}

NetworKit::Cover clean_up(const NetworKit::Graph &graph, const NetworKit::Cover &cover,
                          const std::vector<std::string> &args) {
    // Set parameters
	setUpParameters(args);

    // Import graph
    auto graph_map = read_networkit_graph(graph);
    OslomNetGlobal oslom(graph_map);
	setUpStochastics((int) oslom.stubs());

    // Read input cover
    auto communities = read_networkit_cover(cover);

    // Clean up cover
    auto good_modules = oslom.clean_up(communities, graph.upperNodeIdBound());
	oslom.print_statistics(std::cout, good_modules);

    // Create result cover
    NetworKit::Cover result_cover(graph.upperNodeIdBound());
    result_cover.setUpperBound(cover.upperBound());
    graph.forNodes([&](NetworKit::node u){
        for (auto c : good_modules.memberships[u])
            result_cover.addToSubset(c, u);
    });
    return result_cover;
}

int main(int argc, char** argv) {
    std::vector<std::string> arg_strings;
    for (int i = 1; i < argc; ++i) {
        arg_strings.emplace_back(argv[i]);
        std::cout << argv[i] << std::endl;
    }

    NetworKit::EdgeListReader reader('\t', 0);
    NetworKit::Graph graph = reader.read("/home/armin/Code/graphs/com-amazon.ungraph.txt");
//    NetworKit::Cover C = NetworKit::CoverReader{}.read(
//            "/home/armin/Code/graphs/com-amazon.all.dedup.cmty.txt", graph);

    NetworKit::EgoSplitting ego_algo(graph);
    ego_algo.run();
    auto cover = ego_algo.getCover();

    auto cleaned_cover = clean_up(graph, cover, arg_strings);
    
//    cleaned_cover.forEntries([&](NetworKit::node u, std::set<NetworKit::index> cs){
//       for (auto c : cs)
//           std::cout << u << " in " << c << std::endl;
//    });
}