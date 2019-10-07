from egosplit.benchmarks.evaluation.benchmark_metric import BenchmarkMetric
from egosplit.external import calc_NMI
from networkit.community import CoverF1Similarity


class Time(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		return benchmark.get_time()

	@staticmethod
	def get_name():
		return "Running Time"


class NMI(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		return calc_NMI(benchmark.get_graph(), benchmark.get_cover(),
		                benchmark.get_ground_truth())

	@staticmethod
	def get_name():
		return "NMI"


class F1(BenchmarkMetric):
	"""
	F1-Score of cover->ground-truth. Average of the F1-Score of each detected community.
	"""
	@staticmethod
	def get_value(benchmark):
		return calc_F1(benchmark.get_graph(), benchmark.get_cover(),
		               benchmark.get_ground_truth())

	@staticmethod
	def get_name():
		return "F1-Score"


class F1_rev(BenchmarkMetric):
	"""
	F1-Score of ground-truth->cover. Average of the F1-Score of each ground-truth community.
	"""
	@staticmethod
	def get_value(benchmark):
		return calc_F1(benchmark.get_graph(), benchmark.get_ground_truth(),
		               benchmark.get_cover())

	@staticmethod
	def get_name():
		return "F1-Score (reversed)"


def calc_F1(graph, cover, refCover):
	similarity = CoverF1Similarity(graph, cover, refCover).run()
	return similarity.getUnweightedAverage()