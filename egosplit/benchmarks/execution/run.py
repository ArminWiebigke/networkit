import datetime
import os
from collections import OrderedDict, defaultdict
from copy import copy
from enum import Enum

import egosplit.benchmarks.evaluation.benchmark_metric as bm
import egosplit.benchmarks.evaluation.config
from networkit.stopwatch import clockit
from egosplit.benchmarks.evaluation.metrics import write_results_to_file, add_compact_results, \
	print_compact_results
from egosplit.benchmarks.execution.algo_creation import get_ego_algos, get_other_algos, \
	EgoSplitClusteringAlgorithmsConfig, EgoSplitParameterConfig
from egosplit.benchmarks.evaluation.timings import write_timings
from egosplit.benchmarks.execution.graph_creation import get_graphs, GraphSetsConfig
from egosplit.benchmarks.evaluation.ego_net_partition import analyse_ego_net_partitions
from egosplit.benchmarks.evaluation.stream_to_gephi import stream_partition
from egosplit.benchmarks.evaluation.cover_analysis import analyse_cover
from egosplit.benchmarks.execution.cleanup_functions import cleanup_test
from egosplit.benchmarks.data_structures.cover_benchmark import CoverBenchmark
from egosplit.benchmarks.execution.cleanup import CleanUp, CleanUpConfig
from egosplit.benchmarks.data_structures.benchmark_set import BenchmarkSet


class Evaluation(Enum):
	METRICS = 1
	TIMINGS = 2
	COVER = 3
	EGO_NETS = 4
	STREAM_TO_GEPHI = 5
	CLEANUP = 6


def run_benchmark(benchmark_config: BenchmarkSet, iteration):
	""" Run benchmarks given by a config """
	# setLogLevel('INFO')
	append_results = False
	evaluations = [
		Evaluation.METRICS,
		Evaluation.TIMINGS,
		Evaluation.COVER,
	]
	stream_to_gephi = benchmark_config.stream_to_gephi
	if stream_to_gephi:
		evaluations.append(Evaluation.STREAM_TO_GEPHI)

	store_ego_nets = benchmark_config.store_ego_nets
	if stream_to_gephi:
		store_ego_nets = True
	if store_ego_nets:
		evaluations.append(Evaluation.EGO_NETS)

	result_subfolder = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
	result_dir = '{}/{}-{}/{}/'.format(get_result_dir(), result_subfolder,
	                                   benchmark_config.result_dir, iteration)
	os.makedirs(result_dir)
	result_summary = OrderedDict()

	print('Creating Graphs...')
	graphs = get_graphs(benchmark_config[GraphSetsConfig])

	ego_algos = get_ego_algos(benchmark_config[EgoSplitClusteringAlgorithmsConfig],
	                          benchmark_config[EgoSplitParameterConfig],
	                          benchmark_config[CleanUpConfig],
	                          store_ego_nets)
	other_algos = get_other_algos(benchmark_config.other_algos)
	algos = ego_algos + other_algos

	if stream_to_gephi and len(graphs) * len(algos) > 8:
		raise RuntimeError('Too many runs to stream!')

	print('Starting benchmarks...')
	benchmarks = create_benchmarks(graphs, algos)
	write_scores_per_egonet = benchmark_config.score_per_egonet
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
		egosplit.benchmarks.evaluation.config.Time,
		egosplit.benchmarks.evaluation.config.NMI,
		# bm.F1,
		# bm.F1_rev,
	]
	compact_metric_names = [m.get_name() for m in compact_metrics]
	print_compact_results(summary, compact_metric_names, result_dir)


@clockit
def evaluate_result(graphs, benchmarks, evaluations, append, summary, result_dir,
                    write_scores_per_egonet):
	if Evaluation.METRICS in evaluations:
		# Write results
		metrics = [
			egosplit.benchmarks.evaluation.config.Time,
			egosplit.benchmarks.evaluation.config.NMI,
			egosplit.benchmarks.evaluation.config.F1,
			egosplit.benchmarks.evaluation.config.F1_rev,
		]
		metric_names = [m.get_name() for m in metrics]
		metric_results = defaultdict(lambda: dict())
		for benchmark in benchmarks:
			for metric in metrics:
				metric_results[benchmark][metric.get_name()] = metric.get_value(benchmark)
		write_results_to_file(metric_results, result_dir, metric_names, append)
		# Store results compact in table
		add_compact_results(summary, metric_results, metric_names)
	if Evaluation.TIMINGS in evaluations:
		write_timings(benchmarks, result_dir, append)
	if Evaluation.COVER in evaluations:
		analyse_cover(benchmarks, result_dir, append)
	if Evaluation.EGO_NETS in evaluations:
		analyse_ego_net_partitions(benchmarks, result_dir, append, write_scores_per_egonet)
	if Evaluation.STREAM_TO_GEPHI in evaluations:
		stream_partition(graphs, benchmarks)
	if Evaluation.CLEANUP in evaluations:
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
