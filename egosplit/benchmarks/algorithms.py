import tempfile
import os
import subprocess
from networkit.community import EgoSplitting, OLP
from networkit.stopwatch import Timer
from networkit import graphio
from egosplit.external import remove_small_communities, cleanUpOSLOM

home_path = os.path.expanduser("~")
code_path = home_path + "/Code"
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

	def __exit__(self, exc_type, exc_val, exc_tb):
		self.elapsed = self.timer.stop()
		del self.timer


class CoverAlgorithm:
	def __init__(self):
		self.timer = ContextTimer()
		self.cover = None

	def copy(self):
		""" Override this method in subclass if __init__ takes arguments """
		return self.__class__()

	def getTime(self):
		return self.timer.elapsed

	def getCover(self):
		return self.cover

	def run(self, graph):
		raise NotImplementedError("Run method not implemented!")


class EgoSplitAlgorithm(CoverAlgorithm):
	def __init__(self, name, local_partition_algorithm,
	             global_partition_algorithm=None, clean_up=""):
		super().__init__()
		self.local_partition_algorithm = local_partition_algorithm
		self.global_partition_algorithm = global_partition_algorithm
		self.name = name
		self.clean_up = clean_up
		self.executionInfo = None
		self.egoNetPartitions = None

	def copy(self):
		return EgoSplitAlgorithm(self._name, self.local_partition_algorithm,
		                         self.global_partition_algorithm, self.clean_up)

	@property
	def name(self):
		name = 'Ego_' + self._name
		if self.clean_up != "":
			name += "_clean_" + self.clean_up
		return name

	@name.setter
	def name(self, value):
		self._name = value

	def run(self, graph):
		algo = EgoSplitting(graph, self.local_partition_algorithm,
		                    self.global_partition_algorithm)
		with self.timer:
			algo.run()
			cover = algo.getCover()
			if self.clean_up == "OSLOM":
				cover = cleanUpOSLOM(graph, cover)
			self.cover = cover

		self.executionInfo = algo.getExecutionInfo()
		self.egoNetPartitions = algo.getEgoNetPartitions()

	# Output timings
	# timings = a.getTimings()
	# for name in sorted(timings.keys()):
	# 	print((name.decode('ASCII') + ": ").ljust(26) + str(timings[name]/1000))
	# timings_str = ""
	# for name in sorted(timings.keys()):
	# 	timings_str += str(timings[name]/1000000).ljust(21)
	# self.out_file.write(timings_str + '\n')
	# print(timings_str)

	# def getExecutionInfo(self):
	# 	return self.executionInfo

	def getEgoNetPartitions(self):
		return self.egoNetPartitions


class OlpAlgorithm(CoverAlgorithm):
	def __init__(self):
		super().__init__()
		self.name = 'OLP'

	def run(self, graph):
		a = OLP(graph)
		with self.timer:
			a.run()
			self.cover = a.getCover()


class MosesAlgorithm(CoverAlgorithm):
	def __init__(self):
		super().__init__()
		self.name = 'MOSES'

	def run(self, graph):
		with tempfile.TemporaryDirectory() as tempdir:
			graph_filename = os.path.join(tempdir, "network.dat")
			output_filename = os.path.join(tempdir, "output.dat")
			graphio.writeGraph(graph, graph_filename,
			                   fileformat=graphio.Format.EdgeListTabZero)
			with self.timer:
				subprocess.call(
					[code_path + "/MOSES/moses-binary-linux-x86-64", graph_filename,
					 output_filename], stdout=dev_null)
			output_filename = remove_small_communities(output_filename)
			cover = graphio.CoverReader().read(output_filename, graph)
			self.cover = cover


class GceAlgorithm(CoverAlgorithm):
	def __init__(self, alpha=1.5):
		super().__init__()
		self.name = 'GCE'
		self.alpha = alpha

	def copy(self):
		return GceAlgorithm(self.alpha)

	def run(self, graph):
		with tempfile.TemporaryDirectory() as tempdir:
			graph_filename = os.path.join(tempdir, "network.edgelist")
			cover_filename = os.path.join(tempdir, "cover.txt")
			cover_file = open(cover_filename, 'x')
			graphio.writeGraph(graph, graph_filename,
			                   fileformat=graphio.Format.EdgeListSpaceZero)
			with self.timer:
				subprocess.call(
					[code_path + "/GCECommunityFinder/GCECommunityFinderUbuntu910",
					 graph_filename, "4", "0.6", str(self.alpha), ".75"],
					stdout=cover_file)
			cover_file.close()
			self.cover = graphio.CoverReader().read(cover_filename, graph)


class OslomAlgorithm(CoverAlgorithm):
	def __init__(self):
		super().__init__()
		self.name = 'OSLOM'

	def run(self, graph):
		with tempfile.TemporaryDirectory() as tempdir:
			graph_filename = os.path.join(tempdir, 'network.dat')
			graphio.writeGraph(graph, graph_filename,
			                   fileformat=graphio.Format.EdgeListTabZero)
			with self.timer:
				subprocess.call([code_path + '/OSLOM2/oslom_undir', '-r', '4', '-hr', '0',
				                 '-uw', '-f', graph_filename], stdout=dev_null)
			result = graphio.CoverReader().read(
				os.path.join(graph_filename + '_oslo_files', 'tp'), graph)
			self.cover = result
