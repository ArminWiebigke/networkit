from egosplit.external import calc_NMI
from networkit.community import CoverF1Similarity


class BenchmarkMetric:
	"""
	This abstract class represents a metric of a benchmark run. The main method get_value
	takes a benchmark result and returns the metric value.
	"""

	@staticmethod
	def get_value(benchmark):
		raise NotImplementedError

	@staticmethod
	def get_name():
		raise NotImplementedError


class Time(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		return benchmark.get_time()

	@staticmethod
	def get_name():
		return 'Running Time'


class NMI(BenchmarkMetric):
	"""
	Normalized mutual information
	"""
	@staticmethod
	def get_value(benchmark):
		return calc_NMI(benchmark.get_graph(), benchmark.get_cover(),
		                benchmark.get_ground_truth())

	@staticmethod
	def get_name():
		return 'NMI'


class F1(BenchmarkMetric):
	"""
	F1-Score of cover compared to ground-truth. Average of the F1-Score of each detected community.
	"""
	@staticmethod
	def get_value(benchmark):
		return calc_F1(benchmark.get_graph(), benchmark.get_cover(),
		               benchmark.get_ground_truth())

	@staticmethod
	def get_name():
		return 'F1-Score'


class F1_rev(BenchmarkMetric):
	"""
	F1-Score of ground-truth compared to cover. Average of the F1-Score of each ground-truth community.
	"""
	@staticmethod
	def get_value(benchmark):
		return calc_F1(benchmark.get_graph(), benchmark.get_ground_truth(),
		               benchmark.get_cover())

	@staticmethod
	def get_name():
		return 'F1-Score (reversed)'


def calc_F1(graph, cover, refCover):
	similarity = CoverF1Similarity(graph, cover, refCover).run()
	return similarity.getUnweightedAverage()