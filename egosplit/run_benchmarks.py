#!/usr/bin/env python3
from multiprocessing import Pool

from networkit.stopwatch import Timer
from egosplit.benchmark_sets import get_benchmark_configs
from egosplit.benchmarks.run import start_benchmarks


def run_benchmark(b, iteration):
	t = Timer()
	start_benchmarks(b, iteration)
	print("Total time for benchmarks: {}s".format(str(t.stop())[:6]))


if __name__ == '__main__':
	for b in get_benchmark_configs():
		pool = Pool(4)
		iterations = 2
		parameters = [(b, i) for i in range(iterations)]
		print(parameters)
		pool.starmap(run_benchmark, parameters, 1)
