#!/usr/bin/env python3

import os

from egosplit.plot_scripts.read_data import read_data
from egosplit.plot_scripts.config import set_sns_style
from egosplit.plot_scripts.configure_plots import run

dirs = [
	'plots/communities',
	'plots/metrics',
	'plots/ego_partition/ego_metrics',
	'plots/ego_partition/metrics',
]
for dir in dirs:
	if not os.path.exists(dir):
		os.makedirs(dir)

set_sns_style()

print("Reading results...")
data = read_data()

print("Creating plots...")
run(data)