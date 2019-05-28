from egosplit.external import genLFR
import datetime


class BenchGraph:
	"""This class represents an input graph for a benchmark. In addition to the graph
	object, the name of the graph, its creation parameters and the created ground truth
	cover are stored.
	"""
	global_id = 0
	def __init__(self, graph, ground_truth, name, parameters=""):
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
		parameters = ""
		for (key, value) in parameter_dict.items():
			parameters += key + "=" + str(value) + ","
		graph, cover = genLFR(**parameter_dict)
		super(LFRGraph, self).__init__(graph, cover, name, parameters)
