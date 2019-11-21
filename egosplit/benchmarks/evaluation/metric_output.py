from collections import OrderedDict

from egosplit.benchmarks.evaluation.utility import create_line
from egosplit.benchmarks.data_structures.cover_benchmark import CoverBenchmark


# Write the metric results to the output file, one line per benchmark run
def write_results_to_file(benchmark_results, result_dir, metric_names, append):
	if not append:
		print_headers(result_dir, metric_names)
	result_file = open(result_dir + 'metrics.result', 'a')
	for benchmark, metric_results in benchmark_results.items():
		metric_values = [metric_results[metric] for metric in metric_names]
		line = create_line(*benchmark.output_line(), *metric_values)
		result_file.write(line)


# Print output file headers
def print_headers(result_dir, metric_names):
	with open(result_dir + 'metrics.result', 'w') as f:
		f.write(create_line(*CoverBenchmark.output_header(), *metric_names))


# Print a summary of the metric results in a compact form (table)
def print_compact_results(results, metric_names, result_dir):
	output_file = open(result_dir + 'compact.result', 'a')
	print_results_compact(results, metric_names, output_file)


# Add benchmark results to the compact results
def add_compact_results(summary, benchmark_results, metric_names):
	for benchmark, metric_results in benchmark_results.items():
		algo_name = benchmark.get_algo_name()
		if algo_name not in summary:
			summary[algo_name] = OrderedDict()
		algo_dict = summary[algo_name]
		graph_name = benchmark.get_graph_name()
		if graph_name not in algo_dict:
			d = OrderedDict()
			for m in metric_names:
				d[m] = []
			algo_dict[graph_name] = d
		for metric in metric_names:
			algo_dict[graph_name][metric].append(metric_results[metric])
		# Reassign algo_dict so it works with multiprocessing.Manager().dict
		summary[algo_name] = algo_dict


# Print the compact results to the standard output and the output file
def print_results_compact(results, metrics, output_file):
	algos = list(results.keys())
	graphs = list(results[algos[0]].keys())
	result_decimals = 4
	pre_dot = 6
	str_width = max(pre_dot + 1 + result_decimals,
	                *[len(metric) for metric in metrics]) + 1
	str_width_first = max([len(algo) for algo in algos]) + 1

	# Graph header
	graph_header = '\n' + 'Graph '.rjust(str_width_first)
	for graph in graphs:
		graph_header += '| ' + graph.rjust(len(metrics) * str_width - 1) + ' '
	graph_header += '\n' + str().ljust(
		str_width_first + (2 + len(metrics) * str_width) * len(graphs), '-')
	# Metrics header
	graph_header += '\n' + 'Algo'.ljust(str_width_first)
	value_width = str_width - 1
	for _ in graphs:
		graph_header += '| '
		for metric in metrics:
			graph_header += metric.rjust(value_width) + ' '
	graph_header += '\n' + str().ljust(
		str_width_first + (2 + len(metrics) * str_width) * len(graphs), '-')
	output_str = graph_header

	# Results
	for algo in algos:
		algo_results = (algo + ' ').ljust(str_width_first)
		for graph in graphs:
			algo_results += '| '
			for metric in metrics:
				try:
					value = results[algo][graph][metric]
				except KeyError as e:
					value = 'Missing'
					print(algo, graph, metric, e, sep=', ')
				if isinstance(value, str):
					value_string = value
				else:
					assert(isinstance(value, list))
					average = sum(value) / len(value)
					value_string = '{val:.{d}f}'.format(val=average, d=result_decimals)
				algo_results += value_string.rjust(value_width) + ' '
		output_str += '\n' + algo_results
	output_str += '\n'
	print(output_str)
	output_file.write(output_str)


# Convert a float to a string with a fixed number of decimals
def fixed_decimals(val, decimals=3):
	return '{:.{decimals}f}'.format(val, decimals=decimals)
