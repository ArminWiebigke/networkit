
from networkit import graph
from networkit import graphio
from networkit import community
from networkit import structures
import tempfile
import os
import subprocess
import scipy
import random
import igraph
import leidenalg

home_path = os.path.expanduser('~')
code_path = home_path + '/Code'
dev_null = open(os.devnull, 'w')


def clusterBigClam(graph, numClus, minCom = 5, maxCom = 100):
	with tempfile.TemporaryDirectory() as tempdir:
		graphPath = os.path.join(tempdir, 'graph.edgelist.txt')
		outPath = os.path.join(tempdir, 'cmtyvv.txt')
		graphio.writeGraph(graph, graphPath, graphio.Format.EdgeListTabOne)
		subprocess.call([code_path + '/snap/examples/bigclam/bigclam', '-o:' + tempdir + '/', '-i:' + graphPath, '-c:' + str(numClus), '-mc:' + str(minCom), '-xc:' + str(maxCom)])
		return graphio.SNAPEdgeListPartitionReader().read(outPath, [(i + 1, i) for i in range(graph.upperNodeIdBound())], graph)


# http://www.mapequation.org/code.html
def clusterInfomap(G):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, 'network.txt')
		graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListSpaceZero)
		subprocess.call([code_path + '/infomap/Infomap', '-s', str(random.randint(-2**31, 2**31)), '-2', '-z', ('-d' if G.isDirected() else '-u'), '--clu', graph_filename, tempdir])
		result = community.readCommunities(os.path.join(tempdir, 'network.clu'), format='edgelist-s0')
		while result.numberOfElements() < G.upperNodeIdBound():
			result.toSingleton(result.extend())

		return result


def clusterLouvain(G):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, 'network.txt')
		bin_filename = os.path.join(tempdir, 'network.bin')
		partition_filename = os.path.join(tempdir, 'partition.txt')
		partition_file = open(partition_filename, 'x')
		graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListSpaceZero)
		subprocess.call([code_path + '/louvain/convert', '-i', graph_filename, '-o', bin_filename])
		subprocess.call([code_path + '/louvain/community', bin_filename, '-l', '-1'], stdout=partition_file)
		partition_file.close()
		result = community.LouvainPartitionReader().read(partition_filename)
		while result.numberOfElements() < G.upperNodeIdBound():
			result.toSingleton(result.extend())
		return result


# http://oslom.org/software.htm
def clusterOSLOM(G):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, 'network.dat')
		graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListTabZero)
		subprocess.call([code_path + '/OSLOM2/oslom_undir', '-r', '4', '-hr', '0', '-uw', '-f', graph_filename], stdout=dev_null)
		result = graphio.CoverReader().read(os.path.join(graph_filename + '_oslo_files', 'tp'), G)
		return result


def cleanUpOSLOM(G, cover):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, 'network.dat')
		graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListTabZero)
		cover_filename = os.path.join(tempdir, 'cover.dat')
		graphio.CoverWriter().write(cover, cover_filename)
		subprocess.call([code_path + '/OSLOM2/oslom_undir', '-r', '0', '-hr', '0', '-uw',
		                 '-f', graph_filename, '-hint', cover_filename], stdout=dev_null)
		result = graphio.CoverReader().read(os.path.join(graph_filename + '_oslo_files', 'tp'), G)
		return result


# https://sites.google.com/site/aaronmcdaid/downloads
def clusterMOSES(G):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, 'network.dat')
		output_filename = os.path.join(tempdir, 'output.dat')
		graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListTabZero)
		subprocess.call([code_path + '/MOSES/moses-binary-linux-x86-64', graph_filename, output_filename])
		result = graphio.CoverReader().read(output_filename, G)
		return result


# https://sites.google.com/site/greedycliqueexpansion/ (Use Ubuntu 9.10 binary)
def clusterGCE(G, alpha = 1.5):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, 'network.edgelist')
		cover_filename = os.path.join(tempdir, 'cover.txt')
		cover_file = open(cover_filename, 'x')
		graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListSpaceZero)
		subprocess.call([code_path + '/GCECommunityFinder/GCECommunityFinderUbuntu910', graph_filename, '4', '0.6', str(alpha), '.75'], stdout=cover_file)
		cover_file.close()
		C = graphio.CoverReader().read(cover_filename, G)
		# for u in G.nodes():
		# 	if not C.contains(u):
		# 		C.toSingleton(u)
		return C


def writeLeidenPartition(G, partition, filename):
	partition_ids = []
	for i in range(0, partition.n):
		partition_ids.append([])
	for i in range(0, len(partition)):
		for u in partition[i]:
			# if G.hasNode(u):
			partition_ids[u].append(i)

	file = open(filename, 'w')
	for pid in partition_ids:
		if len(pid) > 1:
			raise RuntimeError
		if len(pid) == 1:
			file.write(str(pid[0]))
		file.write("\n")


def convertLeidenPartition(G, la_partition):
	partition = structures.Partition(G.upperNodeIdBound())
	partition.setUpperBound(len(la_partition) + 1)
	for i in range(0, len(la_partition)):
		for u in la_partition[i]:
			if G.hasNode(u):
				if partition.contains(u):
					print("Error")
				partition.addToSubset(i, u)
	return partition


def convertToIgraph(G):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, 'graph.dat')
		graph_writer = graphio.EdgeListWriter(' ', 0)
		graph_writer.write(G, graph_filename)
		graph_i = igraph.Graph.Read_Edgelist(graph_filename)
		return graph_i


# https://github.com/vtraag/leidenalg
def partitionLeiden(G, partition_type_name):
		graph_i = convertToIgraph(G)
		if partition_type_name == "modularity":
			partition_type = leidenalg.ModularityVertexPartition
		elif partition_type_name == "surprise":
			partition_type = leidenalg.SurpriseVertexPartition
		else:
			raise RuntimeError
		la_partition = leidenalg.find_partition(graph_i, partition_type)
		partition = convertLeidenPartition(G, la_partition)
		# TODO: Something with the partition is wrong.

		# Convert to networkit.Partition
		# partition_filename = os.path.join('./partition.dat')
		# writeLeidenPartition(G, partition_i, partition_filename)
		# partition = community.readCommunities(partition_filename, "edgelist-s0")
		# partition = community.PartitionReader().read(partition_filename)
		# partition = readLeidenPartition(partition_filename)
		return partition


# https://sites.google.com/site/andrealancichinetti/files
def genLFR(N=1000, k=25, maxk=50, mu=0.01, t1=2, t2=1, minc=20, maxc=50, on=500, om=3, C=None):
	args = [code_path + '/lfr_graph_generator/benchmark', '-N', N, '-k', k, '-maxk', maxk, '-mu', mu, '-t1', t1, '-t2', t2, '-minc', minc, '-maxc', maxc]
	if on > 0:
		args.extend(['-on', on, '-om', om])
	if C is not None:
		args.extend(['-C', C])

	with tempfile.TemporaryDirectory() as tempdir:
		old_dir = os.getcwd()
		try:
			os.chdir(tempdir)
			with open('time_seed.dat', 'w') as f:
				f.write(str(random.randint(0, 2**31)))
			subprocess.call(map(str, args), stdout=dev_null)
		finally:
			os.chdir(old_dir)

		G = graphio.readGraph(os.path.join(tempdir, 'network.dat'), fileformat=graphio.Format.LFR)
		if on == 0:
			C = community.readCommunities(os.path.join(tempdir, 'community.dat'), format='edgelist-s1')
		else:
			C = graphio.EdgeListCoverReader(1).read(os.path.join(tempdir, 'community.dat'), G)
		return (G, C)


def calc_NMI(graph, cover, ref_cover):
	with tempfile.TemporaryDirectory() as tempdir:
		# https://github.com/aaronmcdaid/Overlapping-NMI
		old_dir = os.getcwd()
		try:
			os.chdir(tempdir)
			with open('cover.dat', 'w') as f:
				communities = cover.upperBound() * ['']
				for u in graph.nodes():
					for s in cover.subsetsOf(u):
						communities[s] += str(u) + ' '
				for com in communities:
					if com is not '':
						f.write(com + '\n')
			with open('ref_cover.dat', 'w') as f:
				communities = ref_cover.upperBound() * ['']
				for u in graph.nodes():
					for s in ref_cover.subsetsOf(u):
						communities[s] += str(u) + ' '
				for com in communities:
					if com is not '':
						f.write(com + '\n')
			try:
				out = subprocess.check_output([code_path + '/Overlapping-NMI/onmi', 'cover.dat', 'ref_cover.dat'])
				out_lines = out.splitlines()
				nmi_line = out_lines[0]
				nmi_val = float(nmi_line.split()[1])
			except subprocess.CalledProcessError:
				nmi_val = 0.0
		finally:
			os.chdir(old_dir)
	if nmi_val > 1.0 or nmi_val < 0.0:
		nmi_val = -1.0
	return nmi_val


def getFacebookData(name, attribute):
	attribute_dict = {
		'student_fac' : 0,
		'gender' : 1,
		'major_index' : 2,
		'second_major' : 3,
		'dorm' : 4,
		'year' : 5,
		'high_school' : 6,
		}

	if attribute not in attribute_dict:
		raise Exception('Attribute {0} not found'.format(attribute))

	fileName = home_path + '/graphs/facebook100/{0}.mat'.format(name)
	matlabObject = scipy.io.loadmat(fileName)
	col = attribute_dict[attribute]
	n = matlabObject['local_info'].shape[0]
	P = structures.Partition(n)
	for u, a in enumerate(matlabObject['local_info'][:,col]):
		a = max(0, a)
		if a >= P.upperBound():
			P.setUpperBound(a+1)
		P.addToSubset(a, u)
	return P


def remove_small_communities(filename):
	original_file = open(filename, 'r')
	cleaned_filename = filename + '.cleaned'
	cleaned_file = open(cleaned_filename, 'w')
	for line in original_file:
		if len(str.split(line)) > 4:
			cleaned_file.write(line)
	return cleaned_filename


def getFacebookGraph(name):
	return graphio.readMat(home_path + '/graphs/facebook100/{0}.mat'.format(name), key='A')


def get_filename(filename, clean):
	if clean:
		filename = remove_small_communities(filename)
	return filename


def getAmazonGraph5000(clean=False):
	g = graphio.readGraph(code_path + '/graphs/com-amazon.ungraph.txt', fileformat=graphio.Format.EdgeListTabZero)
	filename = code_path + '/graphs/com-amazon.top5000.cmty.txt'
	c = graphio.CoverReader().read(get_filename(filename, clean), g)
	return g, c


def getAmazonGraphAll(clean=False):
	g = graphio.readGraph(code_path + '/graphs/com-amazon.ungraph.txt', fileformat=graphio.Format.EdgeListTabZero)
	filename = code_path + '/graphs/com-amazon.all.dedup.cmty.txt'
	c = graphio.CoverReader().read(get_filename(filename, clean), g)
	return g, c


def getDBLPGraph():
	g = graphio.readGraph(code_path + '/graphs/com-dblp.ungraph.txt', fileformat=graphio.Format.EdgeListTabZero)
	c = graphio.CoverReader().read(code_path + '/graphs/com-dblp.top5000.cmty.txt', g)
	return g, c


def getLiveJournalGraph():
	g = graphio.readGraph(code_path + '/graphs/com-lj.ungraph.txt', fileformat=graphio.Format.EdgeListTabZero)
	c = graphio.CoverReader().read(code_path + '/graphs/com-lj.top5000.cmty.txt', g)
	return g, c


def getOrkutGraph():
	g = graphio.readGraph(code_path + '/graphs/com-orkut.ungraph.txt', fileformat=graphio.Format.EdgeListTabZero)
	c = graphio.CoverReader().read(code_path + '/graphs/com-orkut.top5000.cmty.txt', g)
	return g, c