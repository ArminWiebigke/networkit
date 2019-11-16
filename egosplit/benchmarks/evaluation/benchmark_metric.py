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
