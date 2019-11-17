#!/usr/bin/env python3

import os
import sys
import shutil
from multiprocessing import Pool

from egosplit.benchmarks.benchmark_sets import get_benchmark_configs
from egosplit.benchmarks.plot_scripts.bench_config import PlotGraphSetConfig, PlotAlgoSetConfig
from egosplit.benchmarks.plot_scripts.read_data import DataReader
from egosplit.benchmarks.plot_scripts.plot_config import set_sns_style
from egosplit.benchmarks.plot_scripts.create_plots import PlotSetConfig, make_plots
from egosplit.benchmarks.data_structures.benchmark_set import BenchmarkSet

num_args = len(sys.argv)
result_dir = sys.argv[1] if num_args > 1 else '../results/'
plot_dir = sys.argv[2] if num_args > 2 else '../plots/'

if not result_dir[-1] == '/':
	result_dir += '/'


def create_plots(data, output_dir, config: BenchmarkSet):
	assert output_dir[-1] == '/'
	for plot_func in PlotSetConfig.get_plot_functions(config[PlotSetConfig]):
		for graph_set in PlotGraphSetConfig.get_sets(config[PlotGraphSetConfig]):
			for algo_set_name, algo_set in PlotAlgoSetConfig.get_algo_sets(config[PlotAlgoSetConfig]):
				make_plots(config, data, algo_set_name, algo_set, graph_set, output_dir, plot_func)


def run_plots(config: BenchmarkSet):
	set_sns_style()
	filter_subdir = '*{}'.format(config.result_dir)
	data = DataReader(os.path.join(result_dir, filter_subdir))

	plot_sub_dir = plot_dir + config.plot_dir
	print(plot_sub_dir)
	if os.path.exists(plot_sub_dir):
		shutil.rmtree(plot_sub_dir)
	dirs = [
		'metrics',
		'timings',
		'comm_sizes',
		'comm_f1',
		'num_comms',
		'ego_size_metrics',
		'ego_metrics',
	]
	dirs.extend(['tables/' + d for d in dirs])
	for dir in dirs:
		path = os.path.join(plot_sub_dir, dir)
		if not os.path.exists(path):
			os.makedirs(path)

	print('Creating plots...')
	create_plots(data, plot_sub_dir, config)

	for dir in dirs:
		path = os.path.join(plot_sub_dir, dir)
		files = os.listdir(path)
		if not files:
			os.rmdir(path)


if __name__ == '__main__':
	# pool = Pool(1)
	# pool.map(run_plots, get_benchmark_configs(), 1)
	for config in get_benchmark_configs():
		run_plots(config)
