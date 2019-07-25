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
		f.write(create_line(*CoverBenchmark.output_header(), 'Number of Communities'))
	with open(result_dir + 'cover_comm_sizes.result', 'w') as f:
		f.write(create_line(*CoverBenchmark.output_header(), 'Community Size'))
	with open(result_dir + 'cover_node_comms.result', 'w') as f:
		f.write(create_line(*CoverBenchmark.output_header(), 'Number of Communities per Node'))


# Count the number of communities and their sizes
def count_benchmark_cover(result_dir, benchmark):
	count_communities(result_dir, benchmark)


# Count the number of communities and their sizes
def count_communities(result_dir, benchmark):
	cover = benchmark.get_cover()
	# Number of communities
	with open(result_dir + 'cover_num_comms.result', 'a') as f:
		f.write(create_line(*benchmark.output_line(), cover.numberOfSubsets()))

	# Community sizes
	with open(result_dir + 'cover_comm_sizes.result', 'a') as f:
		for u in cover.subsetSizes():
			f.write(create_line(*benchmark.output_line(), log2(u)))

	# Number of Communities per Node
	with open(result_dir + 'cover_node_comms.result', 'a') as f:
		for u in benchmark.get_graph().nodes():
			num_comms = len(cover.subsetsOf(u))
			if num_comms > 0:
				f.write(create_line(*benchmark.output_line(), log2(num_comms)))

		f.write(create_line(*benchmark.output_line(), str(log2(0.5))))  # TODO: Why?
