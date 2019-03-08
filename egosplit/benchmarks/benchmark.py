from collections import OrderedDict

from networkit.community import EgoSplitting, CoverF1Similarity, PLM, PLP, \
	LPPotts, OLP
from egosplit.benchmarks.algorithms import *
from egosplit.benchmarks.cover_analysis import *
from egosplit.benchmarks.cover_benchmark import *
from egosplit.benchmarks.metrics import *
from egosplit.benchmarks.ego_net_partition import *


def start_benchmarks():
	iterations = 5

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
	for om in range(1, 4):
		LFR_graph_args['om_' + str(om)] = {
			'N': 2000, 'k': 18 * om, 'maxk': 120, 'minc': 60, 'maxc': 100,
			't1': 2, 't2': 2, 'mu': 0.2, 'on': 2000, 'om': om}
	# for mu_factor in range(0, 56, 5):
	# 	LFR_graph_args['mu_' + str(mu_factor).rjust(2,'0')] = {
	# 		'N': 1000, 'k': 20, 'maxk': 50, 'minc': 10, 'maxc': 50,
	# 		't1': 2, 't2': 1, 'mu': 0.01 * mu_factor, 'on': 1000, 'om': 2}
	create_LFR_graphs(graphs, LFR_graph_args, iterations)

	# ************************************************************************************
	# *                         Benchmark algorithms                                     *
	# ************************************************************************************
	algos = []
	partition_algos = OrderedDict()
	# partition_algos['PLP'] = [lambda g: PLP(g, 1, 20).run().getPartition()]
	partition_algos['PLM'] = [lambda g: PLM(g, False, 1.0, "none").run().getPartition()]
	# partition_algos['LPPotts'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition()]
	# partition_algos['LPPotts_par'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition(),
	# 								  lambda g: LPPotts(g, 0, 1, 20, True).run().getPartition()]
	# partition_algos['Infomap'] = [lambda g: clusterInfomap(g)]
	# partition_algos['Surprise'] = [lambda g: partitionLeiden(g, "surprise")]
	# partition_algos['Surprise_PLM'] = [lambda g: partitionLeiden(g, "surprise"),
	#                                    lambda g: PLM(g, False, 1.0, "none").run().getPartition()]
	# partition_algos['Surprise_Infomap'] = [
	# 	lambda g: partitionLeiden(g, "surprise"),
	# 	lambda g: clusterInfomap(g)]
	algos += create_egosplit_algorithms(partition_algos)
	algos += create_egosplit_algorithms(partition_algos, clean_up="OSLOM")

	algos.append(OlpAlgorithm())
	algos.append(GceAlgorithm())
	# algos.append(MosesAlgorithm())
	# algos.append(OslomAlgorithm())

	benchmarks = create_benchmarks(graphs, algos)

	# ************************************************************************************
	# *                     Run benchmarks and get results                               *
	# ************************************************************************************
	# ego_file = open(result_dir + 'ego_timings.result', open_mode)
	# count_ground_truth_communities(graphs, result_dir)

	run_benchmarks(benchmarks)

	result_dir = "../results/"
	append = False
	write_results_to_file(benchmarks, result_dir, append)
	print_benchmarks_compact(benchmarks)
	analyse_ego_net_partitions(benchmarks, result_dir, append)
	analyse_cover(graphs, benchmarks, result_dir, append)


def create_LFR_graphs(graphs, graphArgs, iterations):
	for argsName in graphArgs.keys():
		args = graphArgs[argsName]
		graphName = 'LFR_' + argsName
		for i in range(iterations):
			graph_wrapper = LFRGraph(graphName, args)
			graphs.append(graph_wrapper)
	return graphs


def create_egosplit_algorithms(partition_algos, clean_up=""):
	algos = []
	for p_name in partition_algos:
		algos.append(EgoSplitAlgorithm(p_name, *partition_algos[p_name],
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
