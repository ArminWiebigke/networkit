import datetime
import os
from collections import OrderedDict, defaultdict
from copy import copy

import egosplit.benchmarks.evaluation.benchmark_metric as bm
from algo_creation import get_ego_algos, get_other_algos
from evaluation.timings import write_timings
from graph_creation import get_graphs

from networkit.stopwatch import clockit
from networkit import setLogLevel
from .evaluation.metrics import write_results_to_file, add_compact_results, \
	print_compact_results
from .evaluation.ego_net_partition import analyse_ego_net_partitions
from .evaluation.stream_to_gephi import stream_partition
from .evaluation.cover_analysis import analyse_cover
from egosplit.benchmarks.cleanup import cleanup_test
from .cover_benchmark import CoverBenchmark
from .complete_cleanup import CleanUp


def start_benchmarks(benchmark_config):
	# setLogLevel('INFO')
	iterations = 1
	append_results = False
	evaluations = [
		'metrics',
		'timings',
		'cover',
		'ego_nets',
		# 'stream_to_gephi',
		# 'cleanup',
	]
	stream_to_gephi = 'stream_to_gephi' in evaluations
	store_ego_nets = 'ego_nets' in evaluations and benchmark_config.get('store_ego_nets', True)
	if stream_to_gephi:
		store_ego_nets = True

	result_subfolder = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
	result_dir = '{}/{}-{}/'.format(get_result_dir(), result_subfolder, benchmark_config['result_dir'])
	os.makedirs(result_dir)
	result_summary = OrderedDict()

	for iteration in range(iterations):
		print('Creating Graphs...')
		graphs = get_graphs(benchmark_config['graph_sets'], 1)

		ego_algos = get_ego_algos(benchmark_config['ego_part_algos'], benchmark_config['ego_params'],
		                          store_ego_nets)
		other_algos = get_other_algos(benchmark_config.get('other_algos', None))
		algos = ego_algos + other_algos

		if stream_to_gephi and len(graphs) * len(algos) > 8:
			raise RuntimeError('Too many runs to stream!')

		print('Starting benchmarks...')
		benchmarks = create_benchmarks(graphs, algos)
		write_scores_per_egonet = benchmark_config['score_per_egonet']
		if stream_to_gephi:
			run_benchmarks(benchmarks)
			evaluate_result(graphs, benchmarks, evaluations, append_results,
			                result_summary, result_dir, write_scores_per_egonet)
		else:
			for benchmark in benchmarks:
				run_benchmarks([benchmark])
				evaluate_result(graphs, [benchmark], evaluations, append_results,
				                result_summary, result_dir, write_scores_per_egonet)
				benchmark.clear()
				append_results = True

	print_result_summary(result_summary, get_result_dir())


def get_result_dir():
	return './results/'


def print_result_summary(summary, result_dir):
	if not summary:
		return
	compact_metrics = [
		bm.Time,
		bm.NMI,
		# bm.F1,
		# bm.F1_rev,
	]
	compact_metric_names = [m.get_name() for m in compact_metrics]
	print_compact_results(summary, compact_metric_names, result_dir)


@clockit
def evaluate_result(graphs, benchmarks, evaluations, append, summary, result_dir,
                    write_scores_per_egonet):
	if 'metrics' in evaluations:
		# Write results
		metrics = [
			bm.Time,
			bm.NMI,
			bm.F1,
			bm.F1_rev,
		]
		metric_names = [m.get_name() for m in metrics]
		metric_results = defaultdict(lambda: dict())
		for benchmark in benchmarks:
			for metric in metrics:
				metric_results[benchmark][metric.get_name()] = metric.get_value(benchmark)
		write_results_to_file(metric_results, result_dir, metric_names, append)
		# Store results compact in table
		add_compact_results(summary, metric_results, metric_names)
	if 'timings' in evaluations:
		write_timings(benchmarks, result_dir, append)
	if 'cover' in evaluations:
		analyse_cover(benchmarks, result_dir, append)
	if 'ego_nets' in evaluations:
		analyse_ego_net_partitions(benchmarks, result_dir, append, write_scores_per_egonet)
	if 'stream_to_gephi' in evaluations:
		stream_partition(graphs, benchmarks)
	if 'cleanup' in evaluations:
		cleanup_test(benchmarks, result_dir)


def create_benchmarks(graphs, algos_and_clean):
	benchmarks = []
	for graph in graphs:
		for algo, clean_ups in algos_and_clean:
			algo_copy = copy(algo)
			for clean_up in clean_ups:
				benchmarks.append(CoverBenchmark(algo_copy, CleanUp(clean_up), graph))
	return benchmarks


def run_benchmarks(benchmarks):
	for benchmark in benchmarks:
		benchmark.run()
