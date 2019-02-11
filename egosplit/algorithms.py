from networkit.community import EgoSplitting, OLP
from networkit.stopwatch import Timer
from networkit import graphio
import tempfile
import os
import subprocess

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

	def run(self, graph):
		raise NotImplementedError("Run method not implemented")

class EgoSplitAlgorithm(CoverAlgorithm):
	def __init__(self, out_file, name, localPartitionAlgorithm, globalPartitionAlgorithm = None):
		self.out_file = out_file
		self.algorithm = lambda g: EgoSplitting(g, localPartitionAlgorithm,
												globalPartitionAlgorithm)
		self.name = 'EgoSplitting_(' + name + ')'

	def run(self, graph):
		a = self.algorithm(graph)
		t = Timer()
		a.run()
		t.stop()
		self.time = t.elapsed
		self.cover = a.getCover()
		timings = a.getTimings()
		for name in sorted(timings.keys()):
			print(name.decode('ASCII') + ": " + str(timings[name]/1000))
		timings_str = ""
		for name in sorted(timings.keys()):
			timings_str += str(timings[name]/1000000).ljust(21)
		self.out_file.write(timings_str + '\n')
		print(timings_str)


class OLPAlgorithm(CoverAlgorithm):
	def __init__(self):
		self.name = 'OLP'

	def run(self, graph):
		a = OLP(graph)
		t = Timer()
		a.run()
		t.stop()
		self.time = t.elapsed
		self.cover = a.getCover()


class MOSESAlgorithm(CoverAlgorithm):
	def __init__(self):
		self.name = 'MOSES'

	def run(self, graph):
		with tempfile.TemporaryDirectory() as tempdir:
			graph_filename = os.path.join(tempdir, "network.dat")
			output_filename = os.path.join(tempdir, "output.dat")
			graphio.writeGraph(graph, graph_filename, fileformat=graphio.Format.EdgeListTabZero)
			t = Timer()
			subprocess.call([code_path + "/MOSES/moses-binary-linux-x86-64", graph_filename, output_filename], stdout=dev_null)
			t.stop()
			self.time = t.elapsed
			self.cover = graphio.CoverReader().read(output_filename, graph)


class GCEAlgorithm(CoverAlgorithm):
	def __init__(self, alpha = 1.5):
		self.name = 'GCE'
		self.alpha = alpha

	def run(self, graph):
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

	def run(self, graph):
		with tempfile.TemporaryDirectory() as tempdir:
			graph_filename = os.path.join(tempdir, 'network.dat')
			graphio.writeGraph(graph, graph_filename, fileformat=graphio.Format.EdgeListTabZero)
			t = Timer()
			subprocess.call([code_path + '/OSLOM2/oslom_undir', '-r', '4', '-hr', '0', '-uw', '-f', graph_filename], stdout=dev_null)
			t.stop()
			self.time = t.elapsed
			result = graphio.CoverReader().read(os.path.join(graph_filename + '_oslo_files', 'tp'), graph)
			self.cover = result