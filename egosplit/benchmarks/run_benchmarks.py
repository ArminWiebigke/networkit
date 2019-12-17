#!/usr/bin/env python3
import datetime
from multiprocessing import Pool

from networkit.stopwatch import Timer
from egosplit.benchmarks.benchmark_sets import get_benchmark_configs
from egosplit.benchmarks.execution.run import run_benchmark
from networkit import setLogLevel

if __name__ == '__main__':
	t = Timer()
	timestamp = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
	for benchmark_config in get_benchmark_configs():
		run_benchmark(benchmark_config, timestamp)
	print('Total time for benchmarks: {}s'.format(str(t.stop())[:6]))
