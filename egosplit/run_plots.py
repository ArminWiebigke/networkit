#!/usr/bin/env python3

import os
import sys

from egosplit.plot_scripts.read_data import read_data
from egosplit.plot_scripts.config import set_sns_style
from egosplit.plot_scripts.create_plots import run

num_args = len(sys.argv)
result_dir = sys.argv[1] if num_args > 1 else "results"
plot_dir = sys.argv[2] if num_args > 2 else "plots"

dirs = [
	'communities',
	'metrics',
	'ego_partition/ego_metrics',
	'ego_partition/metrics',
]
for dir in dirs:
	dir = os.path.join(plot_dir, dir)
	if not os.path.exists(dir):
		os.makedirs(dir)

set_sns_style()

print("Reading results...")
data = read_data(result_dir)

print("Creating plots...")
run(data, plot_dir)
