from collections import OrderedDict, defaultdict

from networkit.community import EgoSplitting, CoverF1Similarity, PLM, PLP, \
	LPPotts, OLP
from networkit import gephi
from egosplit.benchmarks.algorithms import *
from egosplit.benchmarks.cover_analysis import *
from egosplit.benchmarks.cover_benchmark import *
from egosplit.benchmarks.metrics import *
from egosplit.benchmarks.ego_net_partition import *


def get_graphs(iterations):
	# ************************************************************************************
	# *                             Input graphs                                         *
	# ************************************************************************************
	graphs = []
	LFR_graph_args = OrderedDict()
	# add_graph(graphs, *getAmazonGraph5000(), "real_Amazon_5000")
	# add_graph(graphs, *getAmazonGraph5000(True), "real_Amazon_5000_no_small")
	# add_graph(graphs, *getAmazonGraphAll(), "real_Amazon_All")
	# add_graph(graphs, *getAmazonGraphAll(True), "real_Amazon_All_no_small")
	# add_graph(graphs, *getDBLPGraph(), "real_DBLP")
	# add_graph(graphs, *getLiveJournalGraph(), 'real_LiveJournal')
	# add_graph(graphs, *getOrkutGraph(), 'real_Orkut')
	# LFR_graph_args['0.01'] = {'N': 1000, 'k': 25, 'maxk': 50, 'mu': 0.01, 'minc': 20,
	# 						'maxc': 50, 'on': 500, 'om': 3}
	# LFR_graph_args['0.1'] = {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.1, 'minc': 5,
	# 					   'maxc': 50, 'on': 100, 'om': 2}
	# LFR_graph_args['0.3'] = {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.3, 'minc': 5,
	# 					   'maxc': 50, 'on': 100, 'om': 2}
	for om in range(5, 6):
		LFR_graph_args['om_' + str(om)] = {
			'N': 2000, 'k': 18 * om, 'maxk': 120, 'minc': 60, 'maxc': 100,
			't1': 2, 't2': 2, 'mu': 0.2, 'on': 2000, 'om': om}
	# for mu_factor in range(0, 51, 5):
	# 	LFR_graph_args['mu_' + str(mu_factor).rjust(2,'0')] = {
	# 		'N': 1000, 'k': 20, 'maxk': 50, 'minc': 10, 'maxc': 50,
	# 		't1': 2, 't2': 1, 'mu': 0.01 * mu_factor, 'on': 1000, 'om': 2}
	create_LFR_graphs(graphs, LFR_graph_args, iterations)
	return graphs


def get_algos():
	# ************************************************************************************
	# *                         Benchmark algorithms                                     *
	# ************************************************************************************
	algos = []
	# algos.append(OslomAlgorithm())
	# algos.append(MosesAlgorithm())
	# algos.append(OlpAlgorithm())
	# algos.append(GceAlgorithm())

	partition_algos = OrderedDict()
	# partition_algos['PLP'] = [lambda g: PLP(g, 1, 20).run().getPartition()]
	partition_algos['PLM'] = [lambda g: PLM(g, False, 1.0, "none").run().getPartition()]
	# partition_algos['LPPotts'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition()]
	# partition_algos['LPPotts_par'] = [
	# 	lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition(),
	# 	lambda g: LPPotts(g, 0, 1, 20, True).run().getPartition()]
	# partition_algos['Infomap'] = [lambda g: clusterInfomap(g)]
	# partition_algos['Surprise'] = [lambda g: partitionLeiden(g, "surprise")]
	# partition_algos['Surprise_PLM'] = [
	# 	lambda g: partitionLeiden(g, "surprise"),
	# 	lambda g: PLM(g, False, 1.0, "none").run().getPartition()]

	ego_parameters = get_ego_parameters()

	algos += create_egosplit_algorithms(partition_algos, ego_parameters)
	# algos += create_egosplit_algorithms(partition_algos, clean_up="OSLOM")
	return algos


def get_ego_parameters():
	ego_parameters = OrderedDict()
	ego_parameters['base'] = {
		"extendEgoNet": "No",
		"removeLowDegreeNeig": "No",
		"edgesBetweenNeigNeig": "No",
		"discardNeigOfNeigEdgesAtFirst": "No",
		"discardNonTriangle": "No",
		"searchNeigOfNeighInDirectedGraph": "No",
		"searchEdgesInDirected": "No",
	}
	ego_parameters_standard = {
		"extendEgoNet": "Yes",
		"removeLowDegreeNeig": "Yes",
		"edgesBetweenNeigNeig": "Yes",
		"discardNeigOfNeigEdgesAtFirst": "Yes",
		"minDegreeCleaning": "1",
		"maxDegreeCleaning": "9999999",
		"addNodesFactor": "0.5",
		"discardNonTriangle": "No",
		"searchNeigOfNeighInDirectedGraph": "No",
		"searchEdgesInDirected": "No",
	}
	ego_parameters['standard'] = {
	}
	ego_parameters['directedNodes'] = {
		"searchNeigOfNeighInDirectedGraph": "Yes",
	}
	ego_parameters['directedEdges'] = {
		"searchEdgesInDirected": "Yes",
	}
	ego_parameters['directedBoth'] = {
		"searchNeigOfNeighInDirectedGraph": "Yes",
		"searchEdgesInDirected": "Yes",
	}
	ego_parameters['discardNonTriangle'] = {
		"discardNonTriangle": "Yes",
	}
	ego_parameters['noNeigOfNeigEdges'] = {
		"edgesBetweenNeigNeig": "No",
	}
	# ego_parameters['noNeigOfNeigEdgesOnlyTriangles'] = {
	# 	"edgesBetweenNeigNeig": "No",
	# 	"discardNonTriangle": "Yes",
	# }
	# ego_parameters['neigOfNeigEdgesAtFirst'] = {
	# 	"discardNeigOfNeigEdgesAtFirst": "No",
	# }

	for key, value in ego_parameters.items():
		ego_parameters[key] = {**ego_parameters_standard, **value}
	return ego_parameters


def start_benchmarks():
	iterations = 3
	graphs = get_graphs(iterations)
	algos = get_algos()
	benchmarks = create_benchmarks(graphs, algos)

	run_benchmarks(benchmarks)

	result_dir = "../results/"
	append = True
	write_results_to_file(benchmarks, result_dir, append)
	print_benchmarks_compact(benchmarks)
	analyse_ego_net_partitions(benchmarks, result_dir, append)
	analyse_cover(graphs, benchmarks, result_dir, append)
	# stream_partition(graphs, benchmarks)


def create_ground_truth_partition(node_id, egonet, ground_truth):
	gt_comms = ground_truth.subsetsOf(node_id)
	partition = dict()
	for u in egonet.nodes():
		comms = ground_truth.subsetsOf(u).intersection(gt_comms)
		if len(comms) == 0:
			color = 1000
		elif len(comms) > 1:
			color = 1001
		else:
			color = list(comms)[0]
		partition[u] = color
	return partition


def mark_direct_neighbors(graph, egonet, node_id):
	neighbors = graph.neighbors(node_id)
	neighbor_map = dict()
	for u in egonet.nodes():
		if u in neighbors:
			neighbor_map[u] = 1
		else:
			neighbor_map[u] = 0
	return neighbor_map


def extend_graph_to_node(graph, u):
	for _ in range(graph.upperNodeIdBound(), u + 1):
		v = graph.addNode()
		graph.removeNode(v)


def stream_partition(graphs, benchmarks):
	workspace_id = 1
	for graph in graphs:
		G = graph.graph
		for _ in range(1):
			node_id = G.randomNode()

			for benchmark in benchmarks:
				if benchmark.graph is not graph \
						or not isinstance(benchmark.algo, EgoSplitAlgorithm):
					continue
				client = gephi.streaming.GephiStreamingClient(
					url='http://localhost:8080/workspace' + str(workspace_id))
				egonet = benchmark.algo.getEgoNet(node_id)
				client.exportGraph(egonet)
				gt_partition = create_ground_truth_partition(node_id, egonet,
				                                             graph.ground_truth)
				client.exportNodeValues(egonet, gt_partition, "ground_truth")
				neighbors = mark_direct_neighbors(graph.graph, egonet, node_id)
				client.exportNodeValues(egonet, neighbors, "neighbors")
				partition = benchmark.algo.getEgoNetPartitions()[node_id]
				client.exportNodeValues(egonet, partition, benchmark.algo.name)
				workspace_id += 1


def create_LFR_graphs(graphs, graphArgs, iterations):
	for argsName in graphArgs.keys():
		args = graphArgs[argsName]
		graphName = 'LFR_' + argsName
		for i in range(iterations):
			graph_wrapper = LFRGraph(graphName, args)
			graphs.append(graph_wrapper)
	return graphs


def create_egosplit_algorithms(partition_algos, ego_parameters, clean_up=""):
	algos = []
	for part_name in partition_algos:
		for para_name, parameters in ego_parameters.items():
			algos.append(EgoSplitAlgorithm(
				part_name + "_" + para_name,
				{key.encode("utf-8"): str(value).encode("utf-8")
				 for (key, value) in parameters.items()},
				*partition_algos[part_name],
				clean_up=clean_up))
	return algos


def create_benchmarks(graphs, algos):
	benchmarks = []
	for algo in algos:
		for graph in graphs:
			benchmarks.append(CoverBenchmark(algo.copy(), graph))
	return benchmarks


def run_benchmarks(benchmarks):
	for benchmark in benchmarks:
		benchmark.run()


def count_ground_truth_communities(graphs, result_dir):
	for graph_wrapper in graphs:
		count_communities(result_dir, "ground_truth", graph_wrapper.name,
		                  graph_wrapper.graph, graph_wrapper.ground_truth)


start_benchmarks()

# from networkit.stopwatch import Timer
# from networkit.sparsification import TriangleEdgeScore
#
#
# def count_triangles(graph):
# 	graph.indexEdges()
# 	triangle = TriangleEdgeScore(graph)
# 	t = Timer()
# 	triangle.run()
# 	t.stop()
# 	print("Time for triangle counting: " + str(t.elapsed))
