from egosplit.external import genLFR
import datetime


class BenchGraph:
	"""This class represents an input graph for a benchmark. In addition to the graph
	object, the name of the graph, its creation parameters and the created ground truth
	cover are stored.
	The actual graph is only created (or read) when it is accessed.
	"""
	global_id = 0

	def __init__(self, name, parameters=None):
		self.name = name
		self.parameters = parameters
		self._graph = None
		self._ground_truth = None
		self.id = '{}({})'.format(datetime.datetime.now().isoformat(), self.global_id)
		self.global_id += 1
		# print('Graph '{}' with {} nodes and {} edges'.format(name, self.graph.numberOfNodes(),
		#                                                      self.graph.numberOfEdges()))

	def create_graph_and_ground_truth(self):
		"""
		Returns graph, ground_truth
		"""
		raise NotImplementedError('Can\'t create graph')

	def set_graph_and_gt(self):
		if not self._graph:
			print('Creating graph {}'.format(self.name))
			self._graph, self._ground_truth = self.create_graph_and_ground_truth()
		assert self._graph

	def clear(self):
		del self._graph
		del self._ground_truth

	@property
	def graph(self):
		self.set_graph_and_gt()
		return self._graph

	@property
	def ground_truth(self):
		self.set_graph_and_gt()
		return self._ground_truth


class ReadGraph(BenchGraph):
	def __init__(self, create_func, name, parameters=None):
		self.create_func = create_func
		super().__init__(name, parameters)

	def create_graph_and_ground_truth(self):
		return self.create_func()


class LFRGraph(BenchGraph):
	"""An input graph created by the LFR graph generator."""
	def __init__(self, name='LFR', parameter_dict=None):
		assert(len(parameter_dict) == len(self.parameter_names()))
		self.lfr_parameters = parameter_dict or {}
		super().__init__(name, parameter_dict)

	def create_graph_and_ground_truth(self):
		return genLFR(**self.lfr_parameters)

	@staticmethod
	def parameter_names():
		return ['N', 'k', 'maxk', 'mu', 't1', 't2', 'minc', 'maxc', 'on', 'om']
