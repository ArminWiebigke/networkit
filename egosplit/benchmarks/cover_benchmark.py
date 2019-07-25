from collections import defaultdict

from graph import LFRGraph


class CoverBenchmark:
	""" This class represents a benchmark instance, consisting of an input graph and the
	algorithm used to create a cover.
	"""

	def __init__(self, algo, clean_up, graph):
		self.algo = algo
		self.clean_up = clean_up
		self.graph = graph
		self.cover = None

	def run(self):
		print("\nGraph: " + self.graph.name + ", Algo: " + self.get_algo_name())

		self.algo.run(self.graph)
		algo_cover = self.algo.get_cover()
		print("Ran algorithm in {:.3f}s".format(self.algo.get_time()))

		self.clean_up.run(self.graph.graph, algo_cover,
		                  self.graph.ground_truth)
		self.cover = self.clean_up.get_cover()
		print("Cleaned up cover in {:.3f}s".format(self.clean_up.get_time()))

	def get_graph(self):
		return self.graph.graph

	def get_graph_name(self):
		return self.graph.name

	def get_graph_id(self):
		return self.graph.id

	def get_ground_truth(self):
		return self.graph.ground_truth

	def get_graph_edges(self):
		return self.graph.graph.numberOfEdges()

	def get_graph_nodes(self):
		return self.graph.graph.numberOfNodes()

	def get_cover(self):
		return self.cover

	def get_time(self):
		return self.algo.get_time() + self.clean_up.get_time()

	def get_timings(self):
		try:
			return self.algo.get_timings()
		except AttributeError:
			return {}

	def get_algo_name(self):
		name = self.algo.name
		if self.clean_up.name:
			name += "_" + self.clean_up.name
		return name

	def clear(self):
		del self.algo
		del self.clean_up
		del self.cover

	@staticmethod
	def output_header():
		return ["Algorithm", "Graph Name", "Graph ID", "Number of Nodes", "Number of Edges"] \
		       + LFRGraph.parameter_names()

	def output_line(self):
		line = [self.get_algo_name(), self.get_graph_name(), self.get_graph_id(),
		        self.get_graph_nodes(), self.get_graph_edges()]
		try:
			line += self.graph.get_lfr_parameters()
		except AttributeError:
			line += ["" for _ in LFRGraph.parameter_names()]
		return line


class MetricCache:
	""" This class is used to cache metric results, so that they are only calculated once.
	"""

	def __init__(self):
		self.cache = defaultdict(lambda: dict())

	def get_metric(self, benchmark, metric):
		if metric not in self.cache[benchmark]:
			self.cache[benchmark][metric] = metric.get_value(benchmark)
		return self.cache[benchmark][metric]

	def clear(self):
		self.cache = defaultdict(lambda: dict())
