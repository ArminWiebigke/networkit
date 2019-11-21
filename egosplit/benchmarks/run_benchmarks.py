#!/usr/bin/env python3
import datetime
from multiprocessing import Pool

from networkit.stopwatch import Timer
from egosplit.benchmarks.benchmark_sets import get_benchmark_configs
from egosplit.benchmarks.execution.run import run_benchmark
from networkit import setLogLevel


def run_benchmark_iteration(benchmark_config, iteration, timestamp):
	run_benchmark(benchmark_config, iteration, timestamp)


if __name__ == '__main__':
	# setLogLevel('INFO')
	t = Timer()
	timestamp = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
	for benchmark_config in get_benchmark_configs():
		iterations = benchmark_config.iterations
		# pool = Pool(1)
		# parameters = [(b, i) for i in range(iterations)]
		# print(parameters)
		# pool.starmap(run_benchmark_iteration, parameters, 1)
		for i in range(iterations):
			run_benchmark_iteration(benchmark_config, i, timestamp)
	print('Total time for benchmarks: {}s'.format(str(t.stop())[:6]))
