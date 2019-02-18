from egosplit.external import *
from networkit.community import EgoSplitting, CoverF1Similarity, PLM, PLP, \
	LPPotts, OLP
from networkit import setLogLevel, setNumberOfThreads
from collections import OrderedDict
from egosplit.algorithms import *
from egosplit.graph import *
from math import log2

def calc_F1(graph, cover, refCover):
	similarity = CoverF1Similarity(graph, cover, refCover).run()
	return similarity.getUnweightedAverage()


def create_LFR_graphs(graphs, graphArgs, iterations):
	for argsName in graphArgs.keys():
		args = graphArgs[argsName]
		graphName = 'LFR_' + argsName
		for i in range(iterations):
			graph_wrapper = LFRGraph(graphName, args)
			count_communities("ground_truth", graphName, graph_wrapper.graph,
							  graph_wrapper.ground_truth)
			graphs.append(graph_wrapper)
	return graphs


def fixed_decimals(val, decimals = 3):
	string = str(val)
	split_strs = string.split('.')
	split_strs[1] = split_strs[1].ljust(decimals, '0')
	return split_strs[0] + '.' + split_strs[1][:decimals]


def print_results_to_file(results, file, print_head):
	algos = list(results.keys())
	graphs = list(results[algos[0]].keys())
	metrics = results[algos[0]][graphs[0]].keys()
	result_decimals = 4
	str_width = max(4 + result_decimals, *[len(metric) for metric in metrics]) + 1
	algo_width = max([len(algo) for algo in algos]) + 2
	graph_width = max([len(graph) for graph in graphs]) + 2

	if print_head:
		header = 'algo'.ljust(algo_width) + 'graph'.ljust(graph_width)
		for metric in metrics:
			header += metric.ljust(str_width)
		file.write(header + '\n')
		print(header)
	for algo in algos:
		for graph in graphs:
			for i in range(len(results[algo][graph]['time'])):
				result_line = algo.ljust(algo_width) + graph.ljust(graph_width)
				for metric in metrics:
					result_line += fixed_decimals(results[algo][graph][metric][i],
												  result_decimals).ljust(str_width)
				file.write(result_line + '\n')
				print(result_line)


def print_results_compact(results):
	algos = list(results.keys())
	graphs = list(results[algos[0]].keys())
	metrics = results[algos[0]][graphs[0]].keys()
	result_decimals = 4
	str_width = max(4 + result_decimals, *[len(metric) for metric in metrics]) + 3
	str_width_first = max([len(algo) for algo in algos]) + 1

	# Header
	graph_header = "\n" + "Graph ".ljust(str_width_first)
	for graph in graphs:
		graph_header += "| " + graph.rjust(len(metrics) * str_width - 3) + ' '
	graph_header += '\n' + str().ljust(str_width_first
									   + len(metrics) * str_width * len(graphs), '-')
	graph_header += '\n' + str().ljust(str_width_first)
	for _ in graphs:
		for metric in metrics:
			graph_header += "| " + metric.rjust(str_width - 3) + ' '
	graph_header += '\n' + str().ljust(
		str_width_first + len(metrics) * str_width * len(graphs), '-')
	print(graph_header)

	# Results
	for algo in algos:
		algo_results = (algo + " ").ljust(str_width_first)
		for graph in graphs:
			for metric in metrics:
				r = results[algo][graph][metric]
				val = sum(r) / len(r)
				algo_results += "| " + str(
					fixed_decimals(val, result_decimals)).rjust(str_width - 3) + ' '
		print(algo_results)


def count_communities(algo_name, graph_name, graph, cover):
	num_comms_file = open('cover_num_comms.txt', 'a')
	comm_sizes_file = open('cover_comm_sizes.txt', 'a')
	node_comms_file = open('cover_node_comms.txt', 'a')
	# Number of communities
	num_comms_file.write(algo_name + ' ' + graph_name + ' '
						 + str(cover.numberOfSubsets()) + '\n')

	# Community sizes
	# max_size = 1000
	# size_cnts = [0] * (max_size + 1)
	# for u in cover.subsetSizes():
	# 	if u <= max_size:
	# 		size_cnts[u] += 1
	# for s in range(max_size + 1):
	# 	cnt = size_cnts[s]
	# 	comm_sizes_file.write(algo_name + ' ' + graph_name + ' ' + str(s) + ' '
	# 						  + str(cnt) + '\n')
	for u in cover.subsetSizes():
		comm_sizes_file.write(algo_name + ' ' + graph_name + ' '
							  + str(log2(u)) + '\n')

	# Number of communities per node
	# max_comms = 100
	# comm_cnts = [0] * (max_comms + 1)
	# for u in range(graph.upperNodeIdBound()):
	# 	num_comms = len(cover.subsetsOf(u))
	# 	comm_cnts[num_comms] += 1
	# for s in range(max_comms + 1):
	# 	cnt = comm_cnts[s]
	# 	node_comms_file.write(algo_name + ' ' + graph_name + ' ' + str(s) + ' '
	# 						  + str(cnt) + '\n')
	node_comms_file.write(algo_name + ' ' + graph_name + ' '
						  + str(log2(0.5)) + '\n')
	for u in range(graph.upperNodeIdBound()):
		num_comms = len(cover.subsetsOf(u))
		if num_comms > 0:
			node_comms_file.write(algo_name + ' ' + graph_name + ' '
								  + str(log2(num_comms)) + '\n')


def print_partition_counts(algo_name, graph_name, partition_counts):
	out_file = open('partition_counts.txt', 'a')
	for count in sorted(partition_counts.keys()):
		if "Components" in count.decode('ASCII'):
			out_algo_name = "components"
		else:
			out_algo_name = algo_name
		out_file.write(out_algo_name + " " + graph_name + " " + count.decode('ASCII') + " "
					   + str(partition_counts[count]) + '\n')


def run_benchmarks(algos, graphs):
	results = OrderedDict()
	for algo in algos:
		algo_name = algo.getName()
		print("\t" + algo_name)
		results[algo_name] = OrderedDict()
		for graph_wrapper in graphs:
			graph_name = graph_wrapper.name
			print("\t\t" + graph_name)
			graph = graph_wrapper.graph
			ground_truth = graph_wrapper.ground_truth
			algo.run(graph)
			cover = algo.getCover()
			f1 = calc_F1(graph, cover, ground_truth)
			f1_rev = calc_F1(graph, ground_truth, cover)
			nmi = calc_NMI(graph, cover, ground_truth)
			t = algo.getTime()
			count_communities(algo_name, graph_name, graph, cover)
			if algo.hasPartitionCounts():
				print_partition_counts(algo_name, graph_name, algo.getPartitionCounts())
			print("\t\t\tF1:" + str(f1) + ", F1_rev:" + str(f1_rev)
				  + ", NMI:" + str(nmi) + ", t:" + str(t))
			if not graph_name in results[algo_name]:
				results[algo_name][graph_name] = OrderedDict()
				r = results[algo_name][graph_name]
				for m in ['F1', 'F1_rev', 'NMI', 'time']:
					r[m] = []
			res_dict = results[algo_name][graph_name]
			res_dict['F1'].append(f1)
			res_dict['F1_rev'].append(f1_rev)
			res_dict['NMI'].append(nmi)
			res_dict['time'].append(t)
	return results


def add_graph(graphs, graph, cover, name):
	count_communities("ground_truth", name, graph, cover)
	graphs.append(BenchGraph(graph, cover, name))


def start_benchmarks():
	# setNumberOfThreads(1)
	iterations = 1
	graphs = []
	LFR_graph_args = OrderedDict()
	algos = []
	partition_algos = OrderedDict()

	append = False
	open_mode = 'w'
	if append:
		open_mode = 'a'
	ego_file = open('ego_timings.txt', open_mode)
	if not append:
		s = 21
		ego_file.write('create EgoNets'.ljust(s) + 'a) assign IDs'.ljust(s)
					   + 'b) find triangles'.ljust(s) + 'c) cluster EgoNet'.ljust(s)
					   + 'd) compact EgoNet'.ljust(s) + 'e) EgoNet subsets'.ljust(s)
					   + 'f) EgoNet subsetCnt'.ljust(s) + 'g) clean up'.ljust(s)
					   + 'split Personas'.ljust(s) + 'connect Personas'.ljust(s)
					   + 'Persona Clustering'.ljust(s) + 'create Cover'.ljust(s)
					   + '\n')

		num_comms_file = open('cover_num_comms.txt', 'w')
		num_comms_file.write('algo graph num_comms\n')
		num_comms_file.close()
		comm_sizes_file = open('cover_comm_sizes.txt', 'w')
		comm_sizes_file.write('algo graph comm_size\n')
		comm_sizes_file.close()
		node_comms_file = open('cover_node_comms.txt', 'w')
		node_comms_file.write('algo graph num_comms\n')
		node_comms_file.close()
		partition_counts_file = open('partition_counts.txt', 'w')
		partition_counts_file.write('algo graph partition_name count\n')
		partition_counts_file.close()

	# ************************************************************************************
	# * Input graphs                                                                     *
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
	for mu_factor in range(5, 86, 10):
		LFR_graph_args['mu_' + str(mu_factor).rjust(2,'0')] = {
			'N': 1000, 'k': 20, 'maxk': 50, 'minc': 10, 'maxc': 50,
			't1': 2, 't2': 1, 'mu': 0.01 * mu_factor, 'on': 1000, 'om': 2}
	create_LFR_graphs(graphs, LFR_graph_args, iterations)

	# ************************************************************************************
	# * Benchmark algorithms                                                             *
	# ************************************************************************************
	# partition_algos['PLP'] = [lambda g: PLP(g, 1, 20).run().getPartition()]
	# partition_algos['PLM'] = [lambda g: PLM(g, False, 1.0, "none").run().getPartition()]
	# partition_algos['LPPotts'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition()]
	# partition_algos['LPPotts_par'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition(),
	# 								  lambda g: LPPotts(g, 0, 1, 20, True).run().getPartition()]
	for partition_algo in partition_algos:
		algos.append(EgoSplitAlgorithm(ego_file, partition_algo, *partition_algos[partition_algo]))
	algos.append(OLPAlgorithm())
	algos.append(GCEAlgorithm())
	algos.append(MOSESAlgorithm())
	# algos.append(OSLOMAlgorithm())

	# ************************************************************************************
	# * Run benchmarks and print results                                                 *
	# ************************************************************************************
	results = run_benchmarks(algos, graphs)

	result_file = open('results.txt', open_mode)
	print_results_to_file(results, result_file, not append)
	print_results_compact(results)

# setLogLevel("INFO")
start_benchmarks()

from networkit.stopwatch import Timer
from networkit.sparsification import TriangleEdgeScore

def countTriangles(graph):
	graph.indexEdges()
	triangle = TriangleEdgeScore(graph)
	t = Timer()
	triangle.run()
	t.stop()
	print("Time for triangle counting: " + str(t.elapsed))
