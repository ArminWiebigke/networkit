from .evaluation.cover_analysis import *
from egosplit.external import *
from networkit.community import CoverF1Similarity


class CoverBenchmark:
	def __init__(self, algo, graph):
		self.algo = algo
		self.graph = graph

	def run(self):
		print("Graph: " + self.graph.name + ", Algo: " + self.algo.name)
		self.algo.run(self.graph.graph)
		print("Time: " + str(self.getTime()) + "\n")

	def getMetric(self, name):
		if name == "time":
			return self.getTime()
		if name == "nmi":
			return self.getNMI()
		if name == "f1":
			return self.getF1()
		if name == "f1_rev":
			return self.getF1_rev()

	def getTime(self):
		return self.algo.getTime()

	def getNMI(self):
		return calc_NMI(self.graph.graph, self.algo.getCover(),
		                self.graph.ground_truth)

	def getF1(self):
		return calc_F1(self.graph.graph, self.algo.getCover(),
		               self.graph.ground_truth)

	def getF1_rev(self):
		return calc_F1(self.graph.graph, self.graph.ground_truth,
		               self.algo.getCover())


def calc_F1(graph, cover, refCover):
	similarity = CoverF1Similarity(graph, cover, refCover).run()
	return similarity.getUnweightedAverage()