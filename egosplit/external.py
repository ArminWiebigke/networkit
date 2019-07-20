import tempfile
import os
import subprocess
from collections import defaultdict

import scipy
import random
import igraph
import leidenalg
import graph_tool

from graph_tool.all import OverlapBlockState

from networkit import graph
from networkit import graphio
from networkit import community
from networkit import structures
from networkit import none


home_path = os.path.expanduser('~')
code_path = home_path + '/Code'
graphs_path = home_path + '/graphs'
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
		seed = random.randint(-2**31, 2**31)
		seed = 7645231
		subprocess.call([code_path + '/infomap/Infomap', '-s', str(seed), '-2', '-z',
		                 ('-d' if G.isDirected() else '-u'), '--clu', graph_filename, tempdir],
		                stdout=dev_null)
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
		subprocess.call([code_path + '/OSLOM2/oslom_undir', '-r', '4', '-hr', '0', '-uw',
		                 '-f', graph_filename], stdout=dev_null)
		result = graphio.CoverReader().read(os.path.join(graph_filename + '_oslo_files', 'tp'), G)
		return result


# def cleanUpOslom(G, cover, merge_bad=False, runs=1, cleanup_strategy='both',
#                  check_minimality=False, simple_cleanup=True, max_extend=2, tolerance=0.1):
# 	with tempfile.TemporaryDirectory() as tempdir:
# 		graph_filename = os.path.join(tempdir, 'network.dat')
# 		graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListTabZero)
# 		cover_filename = os.path.join(tempdir, 'cover.dat')
# 		graphio.CoverWriter().write(cover, cover_filename)
# 		bad_groups_filename = os.path.join(tempdir, 'bad_groups.txt')
# 		params = [code_path + '/OSLOM-clean/oslom_undir', '-r', '0', '-hr', '0', '-uw',
# 		          '-singlet',
# 		          '-f', graph_filename, '-hint', cover_filename,
# 		          '-t', str(tolerance),
# 		          '-cup_runs', str(runs),
# 		          '-cu_strat', cleanup_strategy,
# 		          '-max_extend', str(max_extend),
# 		          '-bad_groups_file', bad_groups_filename
# 		          ]
# 		if merge_bad:
# 			params.append('-merge_bad')
# 		if check_minimality:
# 			params.append('-check_min')
# 		if not simple_cleanup:
# 			params.append('-equiv_cup')
# 		print(params)
# 		subprocess.call(params)
# 		result = graphio.CoverReader().read(os.path.join(graph_filename + '_oslo_files',
# 		                                                 'tp'), G)
# 		bad_groups_file = open(bad_groups_filename, 'r')
# 		bad_groups = []
# 		for group in bad_groups_file:
# 			nodes = group.split(' ')[:-1]  # Last item is the newline character
# 			bad_groups.append([int(u) for u in nodes])
# 		return result, bad_groups


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
def clusterGCE(G, alpha=1.5, min_clique=4):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, 'network.edgelist')
		cover_filename = os.path.join(tempdir, 'cover.txt')
		cover_file = open(cover_filename, 'x')
		graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListSpaceZero)
		subprocess.call([code_path + '/GCECommunityFinder/GCECommunityFinderUbuntu910',
		                 graph_filename, str(min_clique), '0.6', str(alpha), '.75'],
		                stdout=cover_file)
		cover_file.close()
		C = graphio.CoverReader().read(cover_filename, G)
		# for u in G.nodes():
		# 	if not C.contains(u):
		# 		C.toSingleton(u)
		return C


def convertCoverToPartition(G, cover):
	partition = structures.Partition(G.upperNodeIdBound())
	partition.setUpperBound(G.upperNodeIdBound())
	singleton_idx = cover.upperBound()

	for u in G.nodes():
		comms = cover.subsetsOf(u)
		if comms:
			# partition.addToSubset(comms.pop(), u)  # TODO: Use best community
			comm_cnts = defaultdict(lambda: 0)
			for v in G.neighbors(u):
				v_comms = cover.subsetsOf(v).intersection(comms)
				for c in v_comms:
					comm_cnts[c] += 1
			best_comm, _ = max([(val, c) for (c, val) in comm_cnts.items()])
			# print(best_comm)
			partition.addToSubset(best_comm, u)

		else:
			partition.addToSubset(singleton_idx, u)
			singleton_idx += 1

	return partition


# https://github.com/vtraag/leidenalg
def partitionLeiden(G, partition_type_name):
	try:
		graph_i = convertToIgraph(G)
		# checkIgraph(G, graph_i)
		if partition_type_name == 'modularity':
			partition_type = leidenalg.ModularityVertexPartition
		elif partition_type_name == 'surprise':
			partition_type = leidenalg.SurpriseVertexPartition
		elif partition_type_name == 'significance':
			partition_type = leidenalg.SignificanceVertexPartition
		else:
			raise RuntimeError
		la_partition = leidenalg.find_partition(graph_i, partition_type)
		partition = convertLeidenPartition(G, la_partition, graph_i)
		# TODO: Something with the partition is wrong.
	except Exception as e:
		print(e)
		exit(1)
	return partition


def leidenSignificance(G):
	graph_i = convertToIgraph(G)
	optimiser = leidenalg.Optimiser()
	profile = optimiser.resolution_profile(graph_i, leidenalg.CPMVertexPartition,
	                                       resolution_range=(0,1), linear_bisection=False,
	                                       min_diff_resolution=0.01)

	best_p = None
	best_sig = -1
	for p in profile:
		# print(p.resolution_parameter)
		sig = leidenalg.SignificanceVertexPartition.FromPartition(p).quality()
		if sig > best_sig:
			best_sig = sig
			best_p = p
	# print(best_p)
	partition = convertLeidenPartition(G, best_p, graph_i)
	return partition


def convertToIgraph(G):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, 'graph.dat')
		# graph_writer = graphio.EdgeListWriter(' ', 0)
		graph_writer = graphio.GMLGraphWriter()
		graph_writer.write(G, graph_filename)
		graph_i = igraph.Graph.Read_GML(graph_filename)
		return graph_i


def convertLeidenPartition(G, la_partition, graph_i):
	partition = structures.Partition(G.upperNodeIdBound())
	partition.setUpperBound(len(la_partition) + 1)
	for i, nodes in enumerate(la_partition):
		for u in nodes:
			u = int(graph_i.vs[u]['id'])
			partition.addToSubset(i, u)
	return partition


def clusterPeacock(G, part_algorithm):
	with tempfile.TemporaryDirectory() as tempdir:
		old_dir = os.getcwd()
		try:
			os.chdir(tempdir)
			graph_filename = os.path.join(tempdir, 'graph.txt')
			graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListSpaceZero)
			subprocess.call([code_path + '/conga/java -cp conga.jar CONGA',
			                 graph_filename, '–e',
			                 '-peacock 0.1'])
			with open(tempdir + "/split-graph.txt", "r") as f:
				for line in f:
					print(line[:-1])

		except:
			print('Error')
		finally:
			os.chdir(old_dir)


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
		except:
			print('Error')
		finally:
			os.chdir(old_dir)
		G = graphio.readGraph(os.path.join(tempdir, 'network.dat'), fileformat=graphio.Format.LFR)
		# if on == 0:
		# 	C = community.readCommunities(os.path.join(tempdir, 'community.dat'), format='edgelist-s1')
		# else:
		C = graphio.EdgeListCoverReader(1).read(os.path.join(tempdir, 'community.dat'), G)
		return G, C


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
		nmi_val = 0.0
	return nmi_val


def getFacebookGraph(name, clean=False):
	f_graph = graphio.readMat(graphs_path + '/facebook100/{0}.mat'.format(name), key='A')
	cover = structures.Cover(f_graph.upperNodeIdBound())
	attributes = ['student_fac', 'gender', 'major_index', 'second_major', 'dorm', 'year',
	              'high_school']
	id_offset = 0
	for attribute in attributes:
		partition = getFacebookData(name, attribute)
		cover.setUpperBound(cover.upperBound() + partition.upperBound())
		for u in f_graph.nodes():
			p_id = partition.subsetOf(u)
			if p_id != none:
				cover.addToSubset(id_offset + p_id, u)
		id_offset = cover.upperBound()

	if clean:
		from egosplit.benchmarks.cleanup import remove_small_comms
		cover = remove_small_comms(f_graph, cover)
	return f_graph, cover


def getFacebookData(name, attribute):
	attribute_dict = {
		'student_fac': 0,
		'gender': 1,
		'major_index': 2,
		'second_major': 3,
		'dorm': 4,
		'year': 5,
		'high_school': 6,
	}

	if attribute not in attribute_dict:
		raise Exception('Attribute {0} not found'.format(attribute))

	fileName = graphs_path + '/facebook100/{0}.mat'.format(name)
	matlabObject = scipy.io.loadmat(fileName)
	col = attribute_dict[attribute]
	n = matlabObject['local_info'].shape[0]
	P = structures.Partition(n)
	for u, a in enumerate(matlabObject['local_info'][:, col]):
		a = max(0, a)
		if a == 0:
			continue
		if a >= P.upperBound():
			P.setUpperBound(a + 1)
		P.addToSubset(a, u)
	return P


def get_filename(filename, clean):
	if clean:
		filename = remove_small_communities(filename)
	return filename


def remove_small_communities(filename):
	original_file = open(filename, 'r')
	cleaned_filename = filename + '.cleaned'
	cleaned_file = open(cleaned_filename, 'w')
	for line in original_file:
		if len(str.split(line)) > 4:
			cleaned_file.write(line)
	return cleaned_filename


# https://snap.stanford.edu/data/
def getAmazonGraph5000(clean=False):
	g = graphio.readGraph(graphs_path + '/com-amazon.ungraph.txt', fileformat=graphio.Format.EdgeListTabZero)
	filename = graphs_path + '/com-amazon.top5000.cmty.txt'
	c = graphio.CoverReader().read(get_filename(filename, clean), g)
	return g, c


def getAmazonGraphAll(clean=False):
	g = graphio.readGraph(graphs_path + '/com-amazon.ungraph.txt', fileformat=graphio.Format.EdgeListTabZero)
	filename = graphs_path + '/com-amazon.all.dedup.cmty.txt'
	c = graphio.CoverReader().read(get_filename(filename, clean), g)
	return g, c


def getDBLPGraph():
	g = graphio.readGraph(graphs_path + '/com-dblp.ungraph.txt', fileformat=graphio.Format.EdgeListTabZero)
	c = graphio.CoverReader().read(graphs_path + '/com-dblp.top5000.cmty.txt', g)
	return g, c


def getLiveJournalGraph():
	g = graphio.readGraph(graphs_path + '/com-lj.ungraph.txt', fileformat=graphio.Format.EdgeListTabZero)
	c = graphio.CoverReader().read(graphs_path + '/com-lj.top5000.cmty.txt', g)
	return g, c


def getOrkutGraph():
	g = graphio.readGraph(graphs_path + '/com-orkut.ungraph.txt', fileformat=graphio.Format.EdgeListTabZero)
	c = graphio.CoverReader().read(graphs_path + '/com-orkut.top5000.cmty.txt', g)
	return g, c


# https://graph-tool.skewed.de/
def calc_entropy(G, cover, **entropy_args):
	with tempfile.TemporaryDirectory() as tempdir:
		graphPath = os.path.join(tempdir, 'graph.txt')
		graphio.writeGraph(G, graphPath, fileformat=graphio.Format.GraphML)
		gtGraph = graph_tool.load_graph(graphPath, 'graphml')

		edge_property = gtGraph.new_edge_property('vector<int>')
		no_cover_bm = OverlapBlockState(gtGraph, edge_property)
		print('No cover:', no_cover_bm.entropy(**entropy_args))

		edge_property_b = gtGraph.new_edge_property('vector<int>')
		for source, target, edge_id in gtGraph.get_edges():
			set_a = set(cover.subsetsOf(source))
			set_b = set(cover.subsetsOf(target))
			overlap = set_a.intersection(set_b)
			if overlap:
				comm_a = list(overlap)[0]
				comm_b = comm_a
				edge_property_b[gtGraph.edge(source, target)] = [comm_a, comm_b]

		block_state_b = OverlapBlockState(gtGraph, edge_property_b)
		entropy = block_state_b.entropy(**entropy_args)
		print('Cover entropy:', entropy)

		# min_bm = graph_tool.inference.minimize.minimize_blockmodel_dl(gtGraph)
		# print('Minimize:', min_bm.entropy(**entropy_args))
	return entropy
