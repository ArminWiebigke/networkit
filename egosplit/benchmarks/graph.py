from egosplit.external import genLFR


class BenchGraph:
	def __init__(self, graph, ground_truth, name, parameters=""):
		self.graph = graph
		self.ground_truth = ground_truth
		self.name = name
		self.parameters = parameters


class LFRGraph(BenchGraph):
	def __init__(self, name="LFR", parameter_dict=None):
		parameters = ""
		for (key, value) in parameter_dict.items():
			parameters += key + "=" + str(value) + ","
		graph, cover = genLFR(**parameter_dict)
		super(LFRGraph, self).__init__(graph, cover, name, parameters)
