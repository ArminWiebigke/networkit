from egosplit.external import calc_NMI, calc_entropy
from networkit.community import CoverF1Similarity


class BenchmarkMetric:
	"""This class represents a metric of a benchmark run. The main method get_value
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
		return "time"


class NMI(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		return calc_NMI(benchmark.get_graph(), benchmark.get_cover(),
		                benchmark.get_ground_truth())

	@staticmethod
	def get_name():
		return "nmi"


class F1(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		return calc_F1(benchmark.get_graph(), benchmark.get_cover(),
		               benchmark.get_ground_truth())

	@staticmethod
	def get_name():
		return "f1"


class F1_rev(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		return calc_F1(benchmark.get_graph(), benchmark.get_ground_truth(),
		               benchmark.get_cover())

	@staticmethod
	def get_name():
		return "f1_rev"


class Entropy(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		e = calc_entropy(benchmark.get_graph(), benchmark.get_cover())
		base = calc_entropy(benchmark.get_graph(), benchmark.get_ground_truth())
		return (base / e) ** 4

	@staticmethod
	def get_name():
		return "entropy"


class Entropy2(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		e = calc_entropy(benchmark.get_graph(), benchmark.get_cover(),
		                 deg_entropy=False)
		base = calc_entropy(benchmark.get_graph(), benchmark.get_ground_truth(),
		                    deg_entropy=False)
		return (base / e) ** 4

	@staticmethod
	def get_name():
		return "entropy2"


class Entropy3(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		e = calc_entropy(benchmark.get_graph(), benchmark.get_cover(),
		                 deg_entropy=False, degree_dl=False)
		base = calc_entropy(benchmark.get_graph(), benchmark.get_ground_truth(),
		                    deg_entropy=False, degree_dl=False)
		return (base / e) ** 4

	@staticmethod
	def get_name():
		return "entropy3"


class Entropy4(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		e = calc_entropy(benchmark.get_graph(), benchmark.get_cover(),
		                 deg_entropy=False, degree_dl=False, edges_dl=False)
		base = calc_entropy(benchmark.get_graph(), benchmark.get_ground_truth(),
		                    deg_entropy=False, degree_dl=False, edges_dl=False)
		return (base / e) ** 4

	@staticmethod
	def get_name():
		return "entropy4"


# Calculate the F1 score
def calc_F1(graph, cover, refCover):
	similarity = CoverF1Similarity(graph, cover, refCover).run()
	return similarity.getUnweightedAverage()
