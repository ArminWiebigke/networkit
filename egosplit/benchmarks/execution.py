from collections import OrderedDict
from math import log2

from networkit.community import EgoSplitting, CoverF1Similarity, PLM, PLP, \
	LPPotts, OLP
from egosplit.benchmarks.graph import *


def print_headers():
	ego_file = open('ego_timings.txt', 'w')
	s = 21
	ego_file.write('create EgoNets'.ljust(s) + 'a) assign IDs'.ljust(s)
				   + 'b) find triangles'.ljust(s) + 'c) cluster EgoNet'.ljust(s)
				   + 'd) compact EgoNet'.ljust(s) + 'e) EgoNet subsets'.ljust(s)
				   + 'f) EgoNet subsetCnt'.ljust(s) + 'g) clean up'.ljust(s)
				   + 'split Personas'.ljust(s) + 'connect Personas'.ljust(s)
				   + 'Persona Clustering'.ljust(s) + 'create Cover'.ljust(s)
				   + '\n')
	ego_file.close()

	num_comms_file = open('cover_num_comms.txt', 'w')
	num_comms_file.write('algo graph num_comms\n')
	num_comms_file.close()
	comm_sizes_file = open('cover_comm_sizes.txt', 'w')
	comm_sizes_file.write('algo graph comm_size\n')
	comm_sizes_file.close()
	node_comms_file = open('cover_node_comms.txt', 'w')
	node_comms_file.write('algo graph num_comms\n')
	node_comms_file.close()
	partition_counts_file = open('execution_info.txt', 'w')
	partition_counts_file.write('algo graph info_name value\n')
	partition_counts_file.close()


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
	if "e" in string:
		return string
	split_strs = string.split('.')
	if '-' not in split_strs[0]:
		split_strs[0] = " " + split_strs[0]
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


def print_execution_info(algo_name, graph_name, partition_counts):
	out_file = open('execution_info.txt', 'a')
	for info_name in sorted(partition_counts.keys()):
		if "Components" in info_name.decode('ASCII'):
			out_algo_name = "components"
		else:
			out_algo_name = algo_name
		out_file.write(out_algo_name + " " + graph_name + " " + info_name.decode('ASCII') + " "
					   + str(partition_counts[info_name]) + '\n')


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
			algo.run(graph, ground_truth)
			cover = algo.getCover()
			f1 = calc_F1(graph, cover, ground_truth)
			f1_rev = calc_F1(graph, ground_truth, cover)
			nmi = calc_NMI(graph, cover, ground_truth)
			if nmi != -1.0 and not 0.0 <= nmi <= 1.0:
				print(nmi)
				raise RuntimeError
			t = algo.getTime()
			count_communities(algo_name, graph_name, graph, cover)
			if algo.hasExecutionInfo():
				print_execution_info(algo_name, graph_name, algo.getExecutionInfo())
			print("\t\t\tF1:" + str(f1) + ", F1_rev:" + str(f1_rev)
				  + ", NMI:" + str(nmi) + ", t:" + str(t))
			if graph_name not in results[algo_name]:
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