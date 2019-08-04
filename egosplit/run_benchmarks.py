#!/usr/bin/env python3
from multiprocessing import Pool

from networkit.stopwatch import Timer
from egosplit.benchmark_sets import get_benchmark_configs
from egosplit.benchmarks.run import start_benchmarks


def run_benchmark(b):
	t = Timer()
	start_benchmarks(b)
	print("Total time for benchmarks: {}s".format(str(t.stop())[:6]))


if __name__ == '__main__':
	pool = Pool(4)
	pool.map(run_benchmark, get_benchmark_configs(), 1)
