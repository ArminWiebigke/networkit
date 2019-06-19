import tempfile
import os
import subprocess
from networkit.community import EgoSplitting, OLP, SLPA
from networkit.stopwatch import Timer
from networkit import graphio
from .cleanup import trim_comms, merge_overlap_comms, entropy_trim, \
	merge_comms_entropy, merge_gt_comms, remove_comms_entropy, remove_small_comms, \
	add_communities
from .complete_cleanup import clean_up_cover

home_path = os.path.expanduser('~')
code_path = home_path + '/Code'
dev_null = open(os.devnull, 'w')


class ContextTimer(object):
	""" Time code like this:
	t = ContextTimer()
	with t:
		<Code>
	print(t.elapsed)
	"""

	def __init__(self):
		self.elapsed = 0.0
		self.timer = None

	def __enter__(self):
		self.timer = Timer()

	# On context exit, add the elapsed time since context enter
	def __exit__(self, exc_type, exc_val, exc_tb):
		self.elapsed += self.timer.stop()
		del self.timer

	def reset(self):
		self.elapsed = 0.0


class CoverAlgorithm:
	"""
	This class represents an algorithm that takes a graph and runs an overlapping
	community detection, return a cover.
	"""

	def __init__(self, name='CoverAlgorithm', clean_up=''):
		self.name = name
		self.clean_up = clean_up
		self.timer = ContextTimer()
		self.graph = None
		self.ground_truth = None
		self.cover = None

	def copy(self):
		""" Override this method in subclass if __init__ takes arguments """
		print(self.__class__)
		return self.__class__(self.name, self.clean_up)

	def get_time(self):
		return self.timer.elapsed

	def get_cover(self):
		return self.cover

	def run(self, graph):
		self.graph = graph.graph
		self.ground_truth = graph.ground_truth
		self.create_cover()
		t_create = self.timer.elapsed
		print("Created cover in {:.3f}s".format(t_create))
		self.clean_up_cover()
		print("Cleaned up cover in {:.3f}s".format(self.timer.elapsed - t_create))


	def create_cover(self):
		"""
		Override this method in subclasses. Set self.cover to the result of the algorithm.
		"""
		raise NotImplementedError('create_cover method not implemented!')

	def clean_up_cover(self):
		with self.timer:
			self.cover = clean_up_cover(self.graph, self.cover, self.ground_truth,
			                            self.clean_up)


class GroundTruth(CoverAlgorithm):
	def __init__(self, name='Ground_Truth', clean_up=''):
		super().__init__(name, clean_up)

	def create_cover(self):
		self.cover = self.ground_truth


class EgoSplitAlgorithm(CoverAlgorithm):
	def __init__(self, name, parameters, local_partition_algorithm,
	             global_partition_algorithm=None, clean_up=''):
		super().__init__(name, clean_up)
		self.local_partition_algorithm = local_partition_algorithm
		self.global_partition_algorithm = global_partition_algorithm
		self.executionInfo = None
		self.egoNetPartitions = None
		self.egoNets = None
		self.parameters = parameters

	def copy(self):
		return EgoSplitAlgorithm(self.name, self.parameters,
		                         self.local_partition_algorithm,
		                         self.global_partition_algorithm, self.clean_up)

	def create_cover(self):
		algo = EgoSplitting(self.graph, self.local_partition_algorithm,
		                    self.global_partition_algorithm)
		algo.setParameters(self.parameters)
		algo.setGroundTruth(self.ground_truth)

		with self.timer:
			algo.run()
			self.cover = algo.getCover()

		self.executionInfo = algo.getExecutionInfo()
		self.egoNetPartitions = algo.getEgoNetPartitions()
		self.egoNets = algo.getEgoNets()

		# Output timings
		timings = algo.getTimings()
		for name in sorted(timings.keys()):
			print(str(int(timings[name]/1000000)).rjust(7) + '  ' + name.decode('ASCII'))
		timings_str = ''
		for name in sorted(timings.keys()):
			timings_str += str(timings[name]/1000000).ljust(21)
		# self.out_file.write(timings_str + '\n')
		# print(timings_str)

	# def getExecutionInfo(self):
	# 	return self.executionInfo

	def ego_net_partition_of(self, u):
		return self.egoNetPartitions[u]

	def ego_net_of(self, u):
		return self.egoNets[u]


class OlpAlgorithm(CoverAlgorithm):
	def __init__(self, name='OLP', clean_up=''):
		super().__init__(name, clean_up)

	def create_cover(self):
		a = OLP(self.graph)
		with self.timer:
			a.run()
			self.cover = a.getCover()


class MosesAlgorithm(CoverAlgorithm):
	def __init__(self, name='MOSES', clean_up=''):
		super().__init__(name, clean_up)

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


class GceAlgorithm(CoverAlgorithm):
	def __init__(self, name='GCE', clean_up='', alpha=1.5, add_name=''):
		super().__init__(name + add_name, clean_up)
		self.alpha = alpha

	def copy(self):
		return GceAlgorithm(self.name, self.clean_up, self.alpha)

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
	def __init__(self, name='OSLOM', clean_up=''):
		super().__init__(name, clean_up)

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
	def __init__(self, name='SLPA', clean_up='', threshold=0.1, numIterations=100):
		super().__init__(name, clean_up)
		self.threshold = threshold
		self.numIterations = numIterations

	def copy(self):
		return SlpaAlgorithm(self.name, self.clean_up, self.threshold, self.numIterations)

	def create_cover(self):
		with self.timer:
			algo = SLPA(self.graph, self.threshold, self.numIterations)
			algo.run()
			self.cover = algo.getCover()
