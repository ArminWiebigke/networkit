from .cover_analysis import *
from .output import create_line


def write_results_to_file(benchmarks, result_dir, append):
	if not append:
		print_headers(result_dir)
	result_file = open(result_dir + 'metrics.result', 'a')
	for benchmark in benchmarks:
		algo = benchmark.algo
		graph_wrapper = benchmark.graph
		graph = graph_wrapper.graph
		cover = algo.getCover()
		ground_truth = graph_wrapper.ground_truth

		time = algo.getTime()
		f1 = calc_F1(graph, cover, ground_truth)
		f1_rev = calc_F1(graph, ground_truth, cover)
		nmi = calc_NMI(graph, cover, ground_truth)

		line = create_line(algo.name, graph_wrapper.name, time, f1, f1_rev,
		                   nmi)
		result_file.write(line)


def print_headers(result_dir):
	with open(result_dir + 'metrics.result', 'w') as f:
		f.write(create_line('algo', 'graph', 'time', 'f1', 'f1_rev', 'nmi'))
	# with open(result_dir + 'execution_info.result', 'w') as f:
	# 	f.write(create_line('algo', 'graph', 'info_name', 'value'))
	# ego_file = open(result_dir + 'ego_timings.result', 'w')
	# s = 21
	# ego_file.write('create EgoNets'.ljust(s) + 'a) assign IDs'.ljust(s)
	#                + 'b) find triangles'.ljust(s) + 'c) cluster EgoNet'.ljust(s)
	#                + 'd) compact EgoNet'.ljust(s) + 'e) EgoNet subsets'.ljust(s)
	#                + 'f) EgoNet subsetCnt'.ljust(s) + 'g) clean up'.ljust(s)
	#                + 'split Personas'.ljust(s) + 'connect Personas'.ljust(s)
	#                + 'Persona Clustering'.ljust(s) + 'create Cover'.ljust(s)
	#                + '\n')
	# ego_file.close()


def print_benchmarks_compact(benchmarks):
	results = get_compact_results(benchmarks)
	print_results_compact(results)


def get_compact_results(benchmarks):
	results = OrderedDict()
	for benchmark in benchmarks:
		algo = benchmark.algo
		graph = benchmark.graph
		algo_name = algo.name
		if algo_name not in results.keys():
			results[algo_name] = OrderedDict()
		graph_name = graph.name
		if graph_name not in results[algo_name].keys():
			d = OrderedDict()
			d['time'] = []
			d['nmi'] = []
			d['f1'] = []
			d['f1_rev'] = []
			results[algo_name][graph_name] = d

		results[algo_name][graph_name]['time'].append(benchmark.getTime())
		results[algo_name][graph_name]['nmi'].append(benchmark.getNMI())
		results[algo_name][graph_name]['f1'].append(benchmark.getF1())
		results[algo_name][graph_name]['f1_rev'].append(benchmark.getF1_rev())
	return results


def print_results_compact(results):
	algos = list(results.keys())
	graphs = list(results[algos[0]].keys())
	metrics = results[algos[0]][graphs[0]].keys()
	result_decimals = 4
	pre_dot = 3
	str_width = max(pre_dot + 1 + result_decimals,
	                *[len(metric) for metric in metrics]) + 2
	str_width_first = max([len(algo) for algo in algos]) + 1

	# Graph header
	graph_header = "\n" + "Graph ".rjust(str_width_first)
	for graph in graphs:
		graph_header += "| " + graph.rjust(len(metrics) * str_width - 2) + ' '
	graph_header += '\n' + str().ljust(str_width_first
									   + len(metrics) * str_width * len(graphs), '-')
	# Metrics header
	graph_header += '\n' + "Algo".ljust(str_width_first)
	for _ in graphs:
		graph_header += "|"
		for metric in metrics:
			graph_header += metric.rjust(str_width - 1) + ' '
	graph_header += '\n' + str().ljust(
		str_width_first + len(metrics) * str_width * len(graphs), '-')
	print(graph_header)

	# Results
	for algo in algos:
		algo_results = (algo + " ").ljust(str_width_first)
		for graph in graphs:
			algo_results += "|"
			for metric in metrics:
				r = results[algo][graph][metric]
				val = sum(r) / len(r)
				algo_results += str(
					fixed_decimals(val, result_decimals)).rjust(str_width - 1) + ' '
		print(algo_results)


def fixed_decimals(val, decimals=3):
	return str("{:.{decimals}f}").format(val, decimals=decimals)