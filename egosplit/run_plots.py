#!/usr/bin/env python3

import os
import sys
import shutil
import subprocess
from multiprocessing import Pool

from egosplit.benchmark_sets import get_benchmark_configs
from egosplit.plot_scripts.read_data import DataReader
from egosplit.plot_scripts.config import set_sns_style
from egosplit.plot_scripts.create_plots import run

num_args = len(sys.argv)
result_dir = sys.argv[1] if num_args > 1 else "results/"
plot_dir = sys.argv[2] if num_args > 2 else "plots/"

if not result_dir[-1] == "/":
	result_dir += "/"

set_sns_style()


def create_plots(config):
	if config['no_plots']:
		return
	print("Reading results...")
	filter_subdir = "*{}".format(config['result_dir'])
	data = DataReader(os.path.join(result_dir, filter_subdir))

	plot_sub_dir = plot_dir + config['plot_dir']
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

	print("Creating plots...")
	run(data, plot_sub_dir, config)

	for dir in dirs:
		path = os.path.join(plot_sub_dir, dir)
		files = os.listdir(path)
		if not files:
			os.rmdir(path)


if __name__ == '__main__':
	pool = Pool(1)
	pool.map(create_plots, get_benchmark_configs(), 1)
# for config in get_benchmark_configs():
# 	pass
