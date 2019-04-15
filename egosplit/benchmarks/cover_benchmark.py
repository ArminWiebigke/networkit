from .evaluation.cover_analysis import *
from egosplit.external import *
from networkit.community import CoverF1Similarity


class CoverBenchmark:
	def __init__(self, algo, graph):
		self.algo = algo
		self.graph = graph

	def run(self):
		print("Graph: " + self.graph.name + ", Algo: " + self.algo.name)
		self.algo.run_with_wrapper(self.graph)
		print("Time: " + str(self.get_time()) + "\n")

	def get_metric(self, name):
		if name == "time":
			return self.get_time()
		if name == "nmi":
			return self.get_nmi()
		if name == "f1":
			return self.get_f1()
		if name == "f1_rev":
			return self.get_f1_rev()
		if name == "entropy":
			return self.get_entropy()

	def get_time(self):
		return self.algo.get_time()

	def get_nmi(self):
		return calc_NMI(self.graph.graph, self.algo.getCover(),
		                self.graph.ground_truth)

	def get_f1(self):
		return calc_F1(self.graph.graph, self.algo.getCover(),
		               self.graph.ground_truth)

	def get_f1_rev(self):
		return calc_F1(self.graph.graph, self.graph.ground_truth,
		               self.algo.getCover())

	def get_entropy(self):
		return calc_entropy(self.graph.graph, self.algo.getCover())


def calc_F1(graph, cover, refCover):
	similarity = CoverF1Similarity(graph, cover, refCover).run()
	return similarity.getUnweightedAverage()