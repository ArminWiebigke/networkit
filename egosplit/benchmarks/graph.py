from egosplit.external import genLFR
import datetime


class BenchGraph:
	"""This class represents an input graph for a benchmark. In addition to the graph
	object, the name of the graph, its creation parameters and the created ground truth
	cover are stored.
	"""
	global_id = 0

	def __init__(self, graph, ground_truth, name, parameters=None):
		print("Graph '{}' with {} nodes and {} edges".format(name, graph.numberOfNodes(),
		                                                   graph.numberOfEdges()))
		self.graph = graph
		self.ground_truth = ground_truth
		self.name = name
		self.parameters = parameters
		self.id = "{}({})".format(datetime.datetime.now().isoformat(), self.global_id)
		self.global_id += 1


class LFRGraph(BenchGraph):
	"""An input graph created by the LFR graph generator."""
	def __init__(self, name="LFR", parameter_dict=None):
		assert(len(parameter_dict) == len(self.parameter_names()))
		self.lfr_paramters = parameter_dict
		graph, cover = genLFR(**parameter_dict)
		super(LFRGraph, self).__init__(graph, cover, name, parameter_dict)

	def get_lfr_parameters(self):
		return [self.lfr_paramters[p] for p in self.parameter_names()]

	@staticmethod
	def parameter_names():
		return ["N", "k", "maxk", "mu", "t1", "t2", "minc", "maxc", "on", "om"]
