from math import log2

from egosplit.benchmarks.cover_benchmark import *
from egosplit.benchmarks.evaluation.output import create_line


# Analyse the result cover of a benchmark run
def analyse_cover(benchmarks, result_dir, append):
	if not append:
		print_headers(result_dir)

	for benchmark in benchmarks:
		count_benchmark_cover(result_dir, benchmark)


# Print output file headers
def print_headers(result_dir):
	with open(result_dir + 'cover_num_comms.result', 'w') as f:
		f.write(create_line('algo', 'graph', 'num_comms'))
	with open(result_dir + 'cover_comm_sizes.result', 'w') as f:
		f.write(create_line('algo', 'graph', 'comm_size'))
	with open(result_dir + 'cover_node_comms.result', 'w') as f:
		f.write(create_line('algo', 'graph', 'num_comms'))


# Count the number of communities and their sizes
def count_benchmark_cover(result_dir, benchmark):
	count_communities(result_dir, benchmark.get_algo_name(), benchmark.get_graph_name(),
	                  benchmark.get_graph(), benchmark.get_cover())


# Count the number of communities and their sizes
def count_communities(result_dir, algo_name, graph_name, graph, cover):
	# Number of communities
	with open(result_dir + 'cover_num_comms.result', 'a') as f:
		f.write(create_line(algo_name, graph_name, cover.numberOfSubsets()))

	# Community sizes
	with open(result_dir + 'cover_comm_sizes.result', 'a') as f:
		for u in cover.subsetSizes():
			f.write(create_line(algo_name, graph_name, log2(u)))

	# Number of communities per node
	with open(result_dir + 'cover_node_comms.result', 'a') as f:
		f.write(algo_name + ' ' + graph_name + ' ' + str(log2(0.5)) + '\n')
		for u in range(graph.upperNodeIdBound()):
			num_comms = len(cover.subsetsOf(u))
			if num_comms > 0:
				f.write(create_line(algo_name, graph_name, log2(num_comms)))
