from egosplit.external import genLFR


class BenchGraph:
	"""This class represents an input graph for a benchmark. In addition to the graph
	object, the name of the graph, its creation parameters and the created ground truth
	cover are stored.
	"""
	def __init__(self, graph, ground_truth, name, parameters=""):
		self.graph = graph
		self.ground_truth = ground_truth
		self.name = name
		self.parameters = parameters


class LFRGraph(BenchGraph):
	"""An input graph created by the LFR graph generator."""
	def __init__(self, name="LFR", parameter_dict=None):
		parameters = ""
		for (key, value) in parameter_dict.items():
			parameters += key + "=" + str(value) + ","
		graph, cover = genLFR(**parameter_dict)
		super(LFRGraph, self).__init__(graph, cover, name, parameters)
