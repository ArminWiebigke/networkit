import tempfile
import os
import subprocess
from copy import copy

from networkit.community import EgoSplitting, OLP, SLPA
from networkit import graphio
from networkit.graph import Graph
from egosplit.benchmarks.data_structures.context_timer import ContextTimer

home_path = os.path.expanduser('~')
code_path = home_path + '/Code'
dev_null = open(os.devnull, 'w')


class CoverAlgorithm:
	"""
	This class represents an algorithm that takes a graph and runs an overlapping
	community detection, return a cover.
	"""

	def __init__(self, name):
		print('Algorithm', name)
		self.name = name
		self.timer = ContextTimer()
		self.graph = None
		self.ground_truth = None
		self.cover = None
		self.has_run = False

	def __copy__(self):
		""" Override this method in subclass if __init__ takes arguments """
		print(self.__class__)
		return self.__class__(self.name)

	def get_time(self):
		return self.timer.elapsed

	def get_cover(self):
		return self.cover

	def run(self, graph):
		if self.has_run:
			return
		self.has_run = True
		self.graph = graph.graph
		self.ground_truth = graph.ground_truth
		self.create_cover()

	def create_cover(self):
		"""
		Override this method in subclasses. Set self.cover to the result of the algorithm.
		"""
		raise NotImplementedError('create_cover method not implemented!')


class GroundTruth(CoverAlgorithm):
	def __init__(self, name='Ground Truth'):
		super().__init__(name)

	def create_cover(self):
		self.cover = self.ground_truth


class EgoSplitAlgorithm(CoverAlgorithm):
	def __init__(self, name, parameters, local_partition_algorithm,
	             global_partition_algorithm=None):
		super().__init__(name)
		# self.output_parameter['Local Clustering Algorithm'], self.local_partition_algorithm = \
		# 	local_partition_algorithm
		# self.output_parameter['Global Clustering Algorithm'], self.global_partition_algorithm = \
		# 	global_partition_algorithm
		self.local_partition_algorithm = \
			local_partition_algorithm
		self.global_partition_algorithm = \
			global_partition_algorithm
		self.executionInfo = None
		self.egoNetPartitions = None
		self.egoNets = None
		self.parameters = parameters
		self.timings = None

	@staticmethod
	def output_parameter_names():
		return ['Local Clustering Algorithm', 'Global Clustering Algorithm',
		        'Extend EgoNet Strategy',
		        'Extend and Partition Iterations', 'Maximum Extend Factor', 'Edges Score Strategy',
		        'connectPersonasStrat', 'signMerge', 'secondarySigExtRounds',
		        'onlyCheckSignOfMaxCandidates', 'Check Candidates Factor', 'onlyUpdatedCandidates']

	def __copy__(self):
		return EgoSplitAlgorithm(self.name, self.parameters,
		                         self.local_partition_algorithm,
		                         self.global_partition_algorithm)

	def create_cover(self):
		algo = EgoSplitting(self.graph, self.local_partition_algorithm,
		                    self.global_partition_algorithm)
		algo.setParameters({key.encode('utf-8'): str(value).encode('utf-8')
		                    for (key, value) in self.parameters.items()})
		algo.setGroundTruth(self.ground_truth)

		with self.timer:
			algo.run()
			self.cover = algo.getCover()

		self.executionInfo = algo.getExecutionInfo()
		self.egoNetPartitions = algo.getEgoNetPartitions()
		self.egoNets = algo.getEgoNets()

		# Output timings
		timings = algo.getTimings()
		self.timings = {k: v / 1e6 for (k, v) in timings.items()}  # Time in ms
		leading_numbers = 0
		for name in sorted(timings.keys()):
			t_name = name.decode('ASCII')
			# print(t_name)
			lead_new = t_name.find(" ")
			if lead_new < leading_numbers:
				print()
			t_name = t_name[lead_new - 1:]
			print((lead_new - 1) * "  " + str(int(timings[name] / 1000000)).rjust(7)
			      + '   ' + t_name)
			leading_numbers = lead_new
		timings_str = ''
		for name in sorted(timings.keys()):
			timings_str += str(timings[name] / 1000000).ljust(21)

	def get_timings(self):
		return self.timings

	# def getExecutionInfo(self):
	# 	return self.executionInfo

	def ego_net_partition_of(self, u):
		return copy(self.egoNetPartitions[u])

	def ego_net_of(self, u):
		graph = Graph(self.graph.upperNodeIdBound())
		try:
			edges = self.egoNets[u]
		except KeyError:
			return graph
		for edge in edges:
			graph.addEdge(edge['u'], edge['v'], edge['weight'])
		for v in graph.nodes():
			if graph.isIsolated(v):
				graph.removeNode(v)
		graph.removeSelfLoops()
		return graph


class OlpAlgorithm(CoverAlgorithm):
	def __init__(self, name='OLP'):
		super().__init__(name)

	def create_cover(self):
		a = OLP(self.graph)
		with self.timer:
			a.run()
			self.cover = a.getCover()


class MosesAlgorithm(CoverAlgorithm):
	def __init__(self, name='MOSES'):
		super().__init__(name)

	def create_cover(self):
		with tempfile.TemporaryDirectory() as tempdir:
			graph_filename = os.path.join(tempdir, 'network.dat')
			output_filename = os.path.join(tempdir, 'output.dat')
			graphio.writeGraph(self.graph, graph_filename,
			                   fileformat=graphio.Format.EdgeListTabZero)
			with self.timer:
				subprocess.call(
					[code_path + '/MOSES/moses-binary-linux-x86-64', graph_filename,
					 output_filename], stdout=dev_null)
			cover = graphio.CoverReader().read(output_filename, self.graph)
			self.cover = cover


class PeacockAlgorithm(CoverAlgorithm):
	def __init__(self, name='Peacock'):
		super().__init__(name)

	def create_cover(self):
		with tempfile.TemporaryDirectory() as tempdir:
			old_dir = os.getcwd()
			try:
				os.chdir(tempdir)
				graph_filename = os.path.join(tempdir, 'graph.txt')
				graphio.writeGraph(self.graph, graph_filename,
				                   fileformat=graphio.Format.EdgeListSpaceZero)
				params = ['java', '-cp', code_path + '/conga/conga.jar', 'CONGA',
				          graph_filename, '-e',
				          # '-peacock', '0.1'
				          ]
				print(params)
				# subprocess.run(params)
				with open(tempdir + "/split-graph.txt", "r") as f:
					for line in f:
						print(line[:-1])
			except Exception as e:
				# print('Error')
				print(e)
				exit(1)
			finally:
				os.chdir(old_dir)


class GceAlgorithm(CoverAlgorithm):
	def __init__(self, name='GCE', alpha=1.5):
		super().__init__(name)
		self.alpha = alpha

	def __copy__(self):
		return GceAlgorithm(self.name, self.alpha)

	def create_cover(self):
		with tempfile.TemporaryDirectory() as tempdir:
			graph_filename = os.path.join(tempdir, 'network.edgelist')
			cover_filename = os.path.join(tempdir, 'cover.txt')
			cover_file = open(cover_filename, 'x')
			graphio.writeGraph(self.graph, graph_filename,
			                   fileformat=graphio.Format.EdgeListSpaceZero)
			with self.timer:
				subprocess.call(
					[code_path + '/GCECommunityFinder/GCECommunityFinderUbuntu910',
					 graph_filename, '4', '0.6', str(self.alpha), '.75'],
					stdout=cover_file)
			cover_file.close()
			self.cover = graphio.CoverReader().read(cover_filename, self.graph)


class GceNetworkitAlgorithm(CoverAlgorithm):
	def __init__(self, name='GCE-NetworKit', alpha=1.5):
		super().__init__(name)
		self.alpha = alpha

	def __copy__(self):
		return GceAlgorithm(self.name, self.alpha)

	def create_cover(self):
		with tempfile.TemporaryDirectory() as tempdir:
			graph_filename = os.path.join(tempdir, 'network.edgelist')
			cover_filename = os.path.join(tempdir, 'cover.txt')
			cover_file = open(cover_filename, 'x')
			graphio.writeGraph(self.graph, graph_filename,
			                   fileformat=graphio.Format.EdgeListSpaceZero)
			with self.timer:
				subprocess.call(
					[code_path + '/GCECommunityFinder/GCECommunityFinderUbuntu910',
					 graph_filename, '4', '0.6', str(self.alpha), '.75'],
					stdout=cover_file)
			cover_file.close()
			self.cover = graphio.CoverReader().read(cover_filename, self.graph)


class OslomAlgorithm(CoverAlgorithm):
	def __init__(self, name='OSLOM'):
		super().__init__(name)

	def create_cover(self):
		with tempfile.TemporaryDirectory() as tempdir:
			graph_filename = os.path.join(tempdir, 'network.dat')
			graphio.writeGraph(self.graph, graph_filename,
			                   fileformat=graphio.Format.EdgeListTabZero)
			with self.timer:
				subprocess.call([code_path + '/OSLOM2/oslom_undir', '-r', '4', '-hr', '0',
				                 '-uw', '-f', graph_filename], stdout=dev_null)
			result = graphio.CoverReader().read(
				os.path.join(graph_filename + '_oslo_files', 'tp'), self.graph)
			self.cover = result


class SlpaAlgorithm(CoverAlgorithm):
	def __init__(self, name='SLPA', threshold=0.1, numIterations=100):
		super().__init__(name)
		self.threshold = threshold
		self.numIterations = numIterations

	def __copy__(self):
		return SlpaAlgorithm(self.name, self.threshold, self.numIterations)

	def create_cover(self):
		with self.timer:
			algo = SLPA(self.graph, self.threshold, self.numIterations)
			algo.run()
			self.cover = algo.getCover()
