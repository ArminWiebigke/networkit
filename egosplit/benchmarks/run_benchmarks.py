#!/usr/bin/env python3
from multiprocessing import Pool

from networkit.stopwatch import Timer
from egosplit.benchmarks.benchmark_sets import get_benchmark_configs
from egosplit.benchmarks.execution.run import run_benchmark


def run_benchmark_iteration(b, iteration):
	run_benchmark(b, iteration)


if __name__ == '__main__':
	t = Timer()
	for b in get_benchmark_configs():
		pool = Pool(1)
		iterations = 1
		parameters = [(b, i) for i in range(iterations)]
		print(parameters)
		pool.starmap(run_benchmark_iteration, parameters, 1)
	print("Total time for benchmarks: {}s".format(str(t.stop())[:6]))
