from collections import OrderedDict

from egosplit.benchmarks.evaluation.cover_analysis import *
from egosplit.benchmarks.evaluation.output import create_line


# Write the metric results to the output file, one line per benchmark run
def write_results_to_file(benchmark_results, result_dir, metrics, append):
	if not append:
		print_headers(result_dir, [m.get_name() for m in metrics])
	result_file = open(result_dir + 'metrics.result', 'a')
	for result in benchmark_results:
		metric_values = []
		for metric in metrics:
			metric_values.append(result.get_metric(metric))

		line = create_line(result.get_algo_name(), result.get_graph_name(), *metric_values)
		result_file.write(line)


# Print output file headers
def print_headers(result_dir, metrics):
	with open(result_dir + 'metrics.result', 'w') as f:
		f.write(create_line('algo', 'graph', *metrics))
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


# Print a summary of the metric results in a compact form (table)
def print_compact_results(results, result_dir):
	output_file = open(result_dir + "compact.results", "a")
	print_results_compact(results, output_file)


# Add benchmark results to the compact results
def add_compact_results(results, benchmark_results, metrics):
	for result in benchmark_results:
		algo_name = result.get_algo_name()
		if algo_name not in results.keys():
			results[algo_name] = OrderedDict()
		graph_name = result.get_graph_name()
		if graph_name not in results[algo_name].keys():
			d = OrderedDict()
			for m in metrics:
				d[m.get_name()] = []
			results[algo_name][graph_name] = d
		for metric in metrics:
			results[algo_name][graph_name][metric.get_name()].append(result.get_metric(metric))


# Print the compact results to the standard output and the output file
def print_results_compact(results, output_file):
	algos = list(results.keys())
	graphs = list(results[algos[0]].keys())
	metrics = results[algos[0]][graphs[0]].keys()
	result_decimals = 4
	pre_dot = 6
	str_width = max(pre_dot + 1 + result_decimals,
	                *[len(metric) for metric in metrics]) + 1
	str_width_first = max([len(algo) for algo in algos]) + 1

	# Graph header
	graph_header = "\n" + "Graph ".rjust(str_width_first)
	for graph in graphs:
		graph_header += "| " + graph.rjust(len(metrics) * str_width - 1) + ' '
	graph_header += '\n' + str().ljust(
		str_width_first + (2 + len(metrics) * str_width) * len(graphs), '-')
	# Metrics header
	graph_header += '\n' + "Algo".ljust(str_width_first)
	for _ in graphs:
		graph_header += "| "
		for metric in metrics:
			graph_header += metric.rjust(str_width - 1) + ' '
	graph_header += '\n' + str().ljust(
		str_width_first + (2 + len(metrics) * str_width) * len(graphs), '-')
	output_str = graph_header

	# Results
	for algo in algos:
		algo_results = (algo + " ").ljust(str_width_first)
		for graph in graphs:
			algo_results += "| "
			for metric in metrics:
				r = results[algo][graph][metric]
				val = sum(r) / len(r)
				algo_results += str(fixed_decimals(val, result_decimals)
				                    ).rjust(str_width - 1) + ' '
		output_str += "\n" + algo_results
	output_str += "\n"
	print(output_str)
	output_file.write(output_str)


# Convert a float to a string with a fixed number of decimals
def fixed_decimals(val, decimals=3):
	return str("{:.{decimals}f}").format(val, decimals=decimals)