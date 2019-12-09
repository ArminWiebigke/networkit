import datetime
import os
from collections import OrderedDict, defaultdict
from copy import copy
from enum import Enum
from multiprocessing import Manager

from networkit import setLogLevel
import egosplit.benchmarks.evaluation.metrics as bm
from networkit.stopwatch import clockit
from egosplit.benchmarks.evaluation.metric_output import write_results_to_file, add_compact_results, \
	print_compact_results
from egosplit.benchmarks.execution.algo_creation import OtherAlgorithms
from egosplit.benchmarks.execution.egosplit_config import EgoSplitClusteringAlgorithmsConfig, EgoSplitParameterConfig, \
	get_ego_algos
from egosplit.benchmarks.evaluation.timings import write_timings
from egosplit.benchmarks.execution.graph_creation import get_graphs, GraphSetsConfig
from egosplit.benchmarks.evaluation.ego_net_partition import analyse_ego_net_partitions
from egosplit.benchmarks.evaluation.stream_to_gephi import stream_partition
from egosplit.benchmarks.evaluation.cover_evaluation import analyze_cover
from egosplit.benchmarks.execution.cleanup_functions import cleanup_test
from egosplit.benchmarks.data_structures.cover_benchmark import CoverBenchmark
from egosplit.benchmarks.execution.cleanup import CleanUp, CleanUpConfig
from egosplit.benchmarks.data_structures.benchmark_set import BenchmarkSet
from egosplit.benchmarks.data_structures.timelimit import run_with_limited_time


class Evaluation(Enum):
	METRICS = 1
	TIMINGS = 2
	COVER = 3
	EGO_NETS = 4
	CLEANUP = 6


def run_benchmark(benchmark_config: BenchmarkSet, iteration, time_stamp):
	""" Run benchmarks given by a config """
	setLogLevel('INFO')
	append_results = False
	evaluations = [
		Evaluation.METRICS,
		Evaluation.TIMINGS,
		Evaluation.COVER,
	]
	stream_to_gephi = benchmark_config.stream_to_gephi

	store_ego_nets = benchmark_config.store_ego_nets
	if stream_to_gephi:
		store_ego_nets = True
	if store_ego_nets:
		evaluations.append(Evaluation.EGO_NETS)
	calc_f1_per_comm = benchmark_config.calc_f1_per_comm
	write_scores_per_egonet = benchmark_config.score_per_egonet
	time_limit = benchmark_config.time_limit

	result_subfolder = time_stamp
	result_dir = '{}/{}-{}/{}/'.format(get_result_dir(), result_subfolder,
	                                   benchmark_config.result_dir, iteration)
	os.makedirs(result_dir)

	result_summary = Manager().dict()
	# result_summary = OrderedDict()

	print('Creating Graphs...')
	graphs = get_graphs(benchmark_config[GraphSetsConfig])

	ego_algos = get_ego_algos(benchmark_config[EgoSplitClusteringAlgorithmsConfig],
	                          benchmark_config[EgoSplitParameterConfig],
	                          benchmark_config[CleanUpConfig],
	                          store_ego_nets)
	other_algos = OtherAlgorithms.get(benchmark_config[OtherAlgorithms])
	algos = other_algos + ego_algos

	if stream_to_gephi and len(graphs) * len(algos) > 8:
		raise RuntimeError('Too many runs to stream!')

	print('Starting benchmarks...')
	if stream_to_gephi:
		benchmarks = create_benchmarks(graphs, algos)
		run_benchmarks(benchmarks)
		stream_partition(graphs, benchmarks)
	else:
		for graph in graphs:
			graph.set_graph_and_gt()  # Create graph now so it can be copied to the processes
			benchmarks = create_benchmarks([graph], algos)

			for benchmark in benchmarks:
				def run_and_evaluate():
					run_benchmarks([benchmark])
					evaluate_result([benchmark], evaluations, append_results,
					                result_summary, result_dir, write_scores_per_egonet, calc_f1_per_comm)

				finished = run_with_limited_time(run_and_evaluate, (), {}, time_limit)
				if finished:
					append_results = True
				else:
					benchmark_did_not_finish(result_summary, benchmark, result_dir)
				benchmark.clear()
			graph.clear()

	print_result_summary(result_summary, get_result_dir())


def benchmark_did_not_finish(result_summary, benchmark, result_dir):
	algo_name = benchmark.get_algo_name()
	graph_name = benchmark.get_graph_name()
	if algo_name not in result_summary:
		result_summary[benchmark.get_algo_name()] = {graph_name: did_not_finish_metrics()}
	else:
		algo_dict = result_summary[benchmark.get_algo_name()]
		algo_dict[graph_name] = did_not_finish_metrics()
		result_summary[benchmark.get_algo_name()] = algo_dict
	with open(result_dir + "did_not_finish.result", "a") as f:
		f.write("{0}, {1}\n".format(graph_name, algo_name))


def get_result_dir():
	return '../results/'


def did_not_finish_metrics():
	metrics = dict()
	for metric in get_compact_metrics():
		metrics[metric.get_name()] = "DNF"
	return metrics


def get_compact_metrics():
	return [
		bm.Time,
		bm.NMI,
		# bm.F1,
		# bm.F1_rev,
	]


def print_result_summary(summary, result_dir):
	if not summary:
		return
	compact_metrics = get_compact_metrics()
	compact_metric_names = [m.get_name() for m in compact_metrics]
	print_compact_results(summary, compact_metric_names, result_dir)


@clockit
def evaluate_result(benchmarks, evaluations, append, summary, result_dir,
                    write_scores_per_egonet, calc_f1_per_comm):
	if Evaluation.METRICS in evaluations:
		evaluate_metrics(append, benchmarks, result_dir, summary)
	if Evaluation.TIMINGS in evaluations:
		write_timings(benchmarks, result_dir, append)
	if Evaluation.COVER in evaluations:
		analyze_cover(benchmarks, result_dir, calc_f1_per_comm, append)
	if Evaluation.EGO_NETS in evaluations:
		analyse_ego_net_partitions(benchmarks, result_dir, append, write_scores_per_egonet)
	if Evaluation.CLEANUP in evaluations:
		cleanup_test(benchmarks, result_dir)


@clockit
def evaluate_metrics(append, benchmarks, result_dir, summary):
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


def create_benchmarks(graphs, algos_and_clean):
	# benchmarks = Manager().list()
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
