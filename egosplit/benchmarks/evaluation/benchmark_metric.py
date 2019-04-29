from egosplit.external import calc_NMI, calc_entropy
from networkit.community import CoverF1Similarity


class BenchmarkMetric:
	def __init__(self):
		raise NotImplementedError("This class can not be instanced!")

	@staticmethod
	def get_value(benchmark):
		raise NotImplementedError

	@staticmethod
	def get_name():
		raise NotImplementedError


class Time(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		return benchmark.algo.get_time()

	@staticmethod
	def get_name():
		return "time"


class NMI(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		return calc_NMI(benchmark.graph.graph, benchmark.algo.get_cover(),
		                benchmark.graph.ground_truth)

	@staticmethod
	def get_name():
		return "nmi"


class F1(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		return calc_F1(benchmark.graph.graph, benchmark.algo.get_cover(),
		               benchmark.graph.ground_truth)

	@staticmethod
	def get_name():
		return "f1"


class F1_rev(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		return calc_F1(benchmark.graph.graph, benchmark.graph.ground_truth,
		       benchmark.algo.get_cover())

	@staticmethod
	def get_name():
		return "f1_rev"


class Entropy(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		e = calc_entropy(benchmark.graph.graph, benchmark.algo.get_cover())
		base = calc_entropy(benchmark.graph.graph, benchmark.graph.ground_truth)
		return (base / e) ** 4

	@staticmethod
	def get_name():
		return "entropy"


class Entropy2(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		e = calc_entropy(benchmark.graph.graph, benchmark.algo.get_cover(),
		                 deg_entropy=False)
		base = calc_entropy(benchmark.graph.graph, benchmark.graph.ground_truth,
		                    deg_entropy=False)
		return (base / e) ** 4

	@staticmethod
	def get_name():
		return "entropy2"


class Entropy3(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		e = calc_entropy(benchmark.graph.graph, benchmark.algo.get_cover(),
		                 deg_entropy=False, degree_dl=False)
		base = calc_entropy(benchmark.graph.graph, benchmark.graph.ground_truth,
		                    deg_entropy=False, degree_dl=False)
		return (base / e) ** 4

	@staticmethod
	def get_name():
		return "entropy3"


class Entropy4(BenchmarkMetric):
	@staticmethod
	def get_value(benchmark):
		e = calc_entropy(benchmark.graph.graph, benchmark.algo.get_cover(),
		                 deg_entropy=False, degree_dl=False, edges_dl=False)
		base = calc_entropy(benchmark.graph.graph, benchmark.graph.ground_truth,
		                    deg_entropy=False, degree_dl=False, edges_dl=False)
		return (base / e) ** 4

	@staticmethod
	def get_name():
		return "entropy4"


def calc_F1(graph, cover, refCover):
	similarity = CoverF1Similarity(graph, cover, refCover).run()
	return similarity.getUnweightedAverage()
