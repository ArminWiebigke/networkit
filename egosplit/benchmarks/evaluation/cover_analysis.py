from math import log2

from egosplit.benchmarks.data_structures.cover_benchmark import *
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
		f.write(create_line(*CoverBenchmark.output_header(), 'Community Size', 'F1 Score'))
	with open(result_dir + 'cover_node_comms.result', 'w') as f:
		f.write(create_line(*CoverBenchmark.output_header(), 'Number of Communities per Node'))


# Count the number of communities and their sizes
def count_benchmark_cover(result_dir, benchmark):
	count_communities(result_dir, benchmark)


# Count the number of communities and their sizes
def count_communities(result_dir, benchmark):
	cover = benchmark.get_cover()
	ground_truth = benchmark.get_ground_truth()
	comm_map = get_communities(benchmark.get_graph(), cover)
	gt_map = get_communities(benchmark.get_graph(), ground_truth)
	comm_sizes = cover.subsetSizeMap()

	# Number of communities
	with open(result_dir + 'cover_num_comms.result', 'a') as f:
		f.write(create_line(*benchmark.output_line(), cover.numberOfSubsets()))

	# Community sizes
	with open(result_dir + 'cover_comm_sizes.result', 'a') as f:
		for u in cover.getSubsetIds():
			comm = comm_map[u]
			size = comm_sizes[u]
			f1 = f1_score(comm, gt_map)
			f.write(create_line(*benchmark.output_line(), log2(size), f1))

	# Number of Communities per Node
	with open(result_dir + 'cover_node_comms.result', 'a') as f:
		for u in benchmark.get_graph().nodes():
			num_comms = len(cover.subsetsOf(u))
			if num_comms > 0:
				f.write(create_line(*benchmark.output_line(), log2(num_comms)))

		# f.write(create_line(*benchmark.output_line(), str(log2(0.5))))  # TODO: Why?


def get_communities(graph, cover):
	comm_map = defaultdict(lambda: set())
	for u in graph.nodes():
		comms = cover.subsetsOf(u)
		for c in comms:
			comm_map[c].add(u)

	return comm_map


def f1_score(community, ground_truth):
	max_f1 = 0.0
	for gt_comm in ground_truth.values():
		overlap = len(gt_comm.intersection(community))
		if overlap == 0:
			continue
		precision = overlap / len(community)
		recall = overlap / len(gt_comm)
		f1 = 2 * precision * recall / (precision + recall)
		max_f1 = max(max_f1, f1)

	return max_f1

















