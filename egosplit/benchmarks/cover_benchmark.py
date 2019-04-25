from .evaluation.cover_analysis import *
from egosplit.external import *
from networkit.community import CoverF1Similarity


class CoverBenchmark:
	def __init__(self, algo, graph):
		self.algo = algo
		self.graph = graph
		self.time = None
		self.nmi = None
		self.f1 = None
		self.f1_rev = None
		self.entropy = None

	def run(self):
		print("\nGraph: " + self.graph.name + ", Algo: " + self.algo.name)
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
		if not self.time:
			self.time = self.algo.get_time()
		return self.time

	def get_nmi(self):
		if not self.nmi:
			self.nmi = calc_NMI(self.graph.graph, self.algo.get_cover(),
			                    self.graph.ground_truth)
		return self.nmi

	def get_f1(self):
		if not self.f1:
			self.f1 = calc_F1(self.graph.graph, self.algo.get_cover(),
			        self.graph.ground_truth)
		return self.f1

	def get_f1_rev(self):
		if not self.f1_rev:
			self.f1_rev = calc_F1(self.graph.graph, self.graph.ground_truth,
			                      self.algo.get_cover())
		return self.f1_rev

	def get_entropy(self):
		if not self.entropy:
			self.entropy = calc_entropy(self.graph.graph, self.algo.get_cover())
		return self.entropy


def calc_F1(graph, cover, refCover):
	similarity = CoverF1Similarity(graph, cover, refCover).run()
	return similarity.getUnweightedAverage()