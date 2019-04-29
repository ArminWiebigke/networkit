from egosplit.external import *
from networkit.community import CoverF1Similarity


class CoverBenchmark:
	def __init__(self, algo, graph):
		self.algo = algo
		self.graph = graph

	def run(self):
		print("\nGraph: " + self.graph.name + ", Algo: " + self.algo.name)
		self.algo.run_with_wrapper(self.graph)
		print("Time: " + str(self.algo.get_time()) + "\n")


class MetricCache:
	def __init__(self, benchmark):
		self.benchmark = benchmark
		self.cached = {}

	def get_algo_name(self):
		return self.benchmark.algo.name

	def get_graph_name(self):
		return self.benchmark.graph.name

	def get_metric(self, metric):
		if metric not in self.cached:
			self.cached[metric] = metric.get_value(self.benchmark)
		return self.cached[metric]


