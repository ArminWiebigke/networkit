import tempfile
import os
import subprocess
from networkit.community import EgoSplitting, OLP
from networkit.stopwatch import Timer
from networkit import graphio
from egosplit.external import remove_small_communities

home_path = os.path.expanduser("~")
code_path = home_path + "/Code"
dev_null = open(os.devnull, 'w')

class CoverAlgorithm:
	cover = None
	time = 0
	name = None

	def getTime(self):
		return self.time

	def getCover(self):
		return self.cover

	def getName(self):
		return self.name

	def run(self, graph, ground_truth):
		raise NotImplementedError("Run method not implemented")

	def hasExecutionInfo(self):
		return False

class EgoSplitAlgorithm(CoverAlgorithm):
	executionInfo = None

	def __init__(self, out_file, name, localPartitionAlgorithm, globalPartitionAlgorithm = None,
				 groundTruth = None):
		self.out_file = out_file
		self.ground_truth = None
		self.algorithm = lambda g, ground_truth: EgoSplitting(g, localPartitionAlgorithm,
												globalPartitionAlgorithm, ground_truth)
		self.name = 'EgoSplitting_(' + name + ')'

	def run(self, graph, ground_truth):
		a = self.algorithm(graph, ground_truth)
		t = Timer()
		a.run()
		t.stop()
		self.time = t.elapsed
		self.cover = a.getCover()
		# Output timings
		timings = a.getTimings()
		for name in sorted(timings.keys()):
			print(name.decode('ASCII') + ": " + str(timings[name]/1000))
		timings_str = ""
		for name in sorted(timings.keys()):
			timings_str += str(timings[name]/1000000).ljust(21)
		self.out_file.write(timings_str + '\n')
		print(timings_str)
		# Output partition counts
		self.executionInfo = a.getExecutionInfo()

	def hasExecutionInfo(self):
		return True

	def getExecutionInfo(self):
		return self.executionInfo


class OLPAlgorithm(CoverAlgorithm):
	def __init__(self):
		self.name = 'OLP'

	def run(self, graph, ground_truth):
		a = OLP(graph)
		t = Timer()
		a.run()
		t.stop()
		self.time = t.elapsed
		self.cover = a.getCover()


class MOSESAlgorithm(CoverAlgorithm):
	def __init__(self):
		self.name = 'MOSES'

	def run(self, graph, ground_truth):
		with tempfile.TemporaryDirectory() as tempdir:
			graph_filename = os.path.join(tempdir, "network.dat")
			output_filename = os.path.join(tempdir, "output.dat")
			graphio.writeGraph(graph, graph_filename, fileformat=graphio.Format.EdgeListTabZero)
			t = Timer()
			subprocess.call([code_path + "/MOSES/moses-binary-linux-x86-64", graph_filename, output_filename], stdout=dev_null)
			t.stop()
			self.time = t.elapsed
			output_filename = remove_small_communities(output_filename)
			cover = graphio.CoverReader().read(output_filename, graph)
			self.cover = cover


class GCEAlgorithm(CoverAlgorithm):
	def __init__(self, alpha = 1.5):
		self.name = 'GCE'
		self.alpha = alpha

	def run(self, graph, ground_truth):
		with tempfile.TemporaryDirectory() as tempdir:
			graph_filename = os.path.join(tempdir, "network.edgelist")
			cover_filename = os.path.join(tempdir, "cover.txt")
			cover_file = open(cover_filename, 'x')
			graphio.writeGraph(graph, graph_filename, fileformat=graphio.Format.EdgeListSpaceZero)
			t = Timer()
			subprocess.call([code_path + "/GCECommunityFinder/GCECommunityFinderUbuntu910", graph_filename, "4", "0.6", str(self.alpha), ".75"], stdout=cover_file)
			t.stop()
			self.time = t.elapsed
			cover_file.close()
			self.cover = graphio.CoverReader().read(cover_filename, graph)


class OSLOMAlgorithm(CoverAlgorithm):
	def __init__(self):
		self.name = 'OSLOM'

	def run(self, graph, ground_truth):
		with tempfile.TemporaryDirectory() as tempdir:
			graph_filename = os.path.join(tempdir, 'network.dat')
			graphio.writeGraph(graph, graph_filename, fileformat=graphio.Format.EdgeListTabZero)
			t = Timer()
			subprocess.call([code_path + '/OSLOM2/oslom_undir', '-r', '4', '-hr', '0', '-uw', '-f', graph_filename], stdout=dev_null)
			t.stop()
			self.time = t.elapsed
			result = graphio.CoverReader().read(os.path.join(graph_filename + '_oslo_files', 'tp'), graph)
			self.cover = result