from egosplit.external import *
from collections import OrderedDict

from networkit.community import EgoSplitting, CoverF1Similarity, PLM, PLP, \
	LPPotts, OLP
from egosplit.benchmarks.algorithms import *
from egosplit.benchmarks.execution import *


def start_benchmarks():
	iterations = 1
	graphs = []
	LFR_graph_args = OrderedDict()
	algos = []
	partition_algos = OrderedDict()

	result_dir = "results/"
	append = False
	if not append:
		print_headers(result_dir)
	open_mode = 'w'
	if append:
		open_mode = 'a'
	ego_file = open(result_dir + 'ego_timings.result', open_mode)

	# ************************************************************************************
	# *                             Input graphs                                         *
	# ************************************************************************************
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
	for om in range(1, 6):
		LFR_graph_args['om_' + str(om)] = {
			'N': 2000, 'k': 18 * om, 'maxk': 120, 'minc': 60, 'maxc': 100,
			't1': 2, 't2': 2, 'mu': 0.2, 'on': 2000, 'om': om}
	# for mu_factor in range(0, 56, 5):
	# 	LFR_graph_args['mu_' + str(mu_factor).rjust(2,'0')] = {
	# 		'N': 1000, 'k': 20, 'maxk': 50, 'minc': 10, 'maxc': 50,
	# 		't1': 2, 't2': 1, 'mu': 0.01 * mu_factor, 'on': 1000, 'om': 2}
	create_LFR_graphs(result_dir, graphs, LFR_graph_args, iterations)

	# ************************************************************************************
	# *                         Benchmark algorithms                                     *
	# ************************************************************************************
	partition_algos['PLP'] = [lambda g: PLP(g, 1, 20).run().getPartition()]
	partition_algos['PLM'] = [lambda g: PLM(g, False, 1.0, "none").run().getPartition()]
	# partition_algos['LPPotts'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition()]
	# partition_algos['LPPotts_par'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition(),
	# 								  lambda g: LPPotts(g, 0, 1, 20, True).run().getPartition()]
	partition_algos['Infomap'] = [lambda g: clusterInfomap(g)]
	partition_algos['Leiden_surprise'] = [lambda g: partitionLeiden(g, "surprise")]
	for partition_algo in partition_algos:
		algos.append(EgoSplitAlgorithm(ego_file, partition_algo, *partition_algos[partition_algo]))
		# algos.append(EgoSplitAlgorithm(ego_file, partition_algo, *partition_algos[partition_algo], clean_up="OSLOM"))
	# algos.append(OLPAlgorithm())
	# algos.append(GCEAlgorithm())
	# algos.append(MOSESAlgorithm())
	algos.append(OSLOMAlgorithm())

	# ************************************************************************************
	# *                     Run benchmarks and print results                             *
	# ************************************************************************************
	results = run_benchmarks(algos, graphs, result_dir)

	result_file = open(result_dir + 'metrics.result', open_mode)
	print_results_to_file(results, result_file, not append)
	print_results_compact(results)

# setLogLevel("INFO")
start_benchmarks()

from networkit.stopwatch import Timer
from networkit.sparsification import TriangleEdgeScore


def count_triangles(graph):
	graph.indexEdges()
	triangle = TriangleEdgeScore(graph)
	t = Timer()
	triangle.run()
	t.stop()
	print("Time for triangle counting: " + str(t.elapsed))
