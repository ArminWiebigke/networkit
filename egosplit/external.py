
from networkit import graph
from networkit import graphio
from networkit import community
from networkit import structures
import tempfile
import os
import subprocess
import scipy
import random

home_path = os.path.expanduser("~")
code_path = home_path + "/Code"
dev_null = open(os.devnull, 'w')

def clusterBigClam(graph, numClus, minCom = 5, maxCom = 100):
	with tempfile.TemporaryDirectory() as tempdir:
		graphPath = os.path.join(tempdir, "graph.edgelist.txt")
		outPath = os.path.join(tempdir, "cmtyvv.txt")
		graphio.writeGraph(graph, graphPath, graphio.Format.EdgeListTabOne)
		subprocess.call([code_path + "/snap/examples/bigclam/bigclam", "-o:" + tempdir + "/", "-i:" + graphPath, "-c:" + str(numClus), "-mc:" + str(minCom), "-xc:" + str(maxCom)])
		return graphio.SNAPEdgeListPartitionReader().read(outPath, [(i + 1, i) for i in range(graph.upperNodeIdBound())], graph)

def clusterInfomap(G):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, "network.txt")
		graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListSpaceZero)
		subprocess.call([code_path + "/infomap/Infomap", "-s", str(random.randint(-2**31, 2**31)), "-2", "-z", ("-d" if G.isDirected() else "-u"), "--clu", graph_filename, tempdir])
		result = community.readCommunities(os.path.join(tempdir, "network.clu"), format="edgelist-s0")
		while result.numberOfElements() < G.upperNodeIdBound():
			result.toSingleton(result.extend())

		return result

def clusterLouvain(G):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, "network.txt")
		bin_filename = os.path.join(tempdir, "network.bin")
		partition_filename = os.path.join(tempdir, "partition.txt")
		partition_file = open(partition_filename, 'x')
		graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListSpaceZero)
		subprocess.call([code_path + "/louvain/convert", "-i", graph_filename, "-o", bin_filename])
		subprocess.call([code_path + "/louvain/community", bin_filename, "-l", "-1"], stdout=partition_file)
		partition_file.close()
		result = community.LouvainPartitionReader().read(partition_filename)
		while result.numberOfElements() < G.upperNodeIdBound():
			result.toSingleton(result.extend())
		return result

def clusterOSLOM(G):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, "network.dat")
		graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListTabZero)
		subprocess.call([code_path + "/OSLOM2/oslom_undir", "-r", "4", "-hr", "0", "-uw", "-f", graph_filename], stdout=dev_null)
		result = graphio.CoverReader().read(os.path.join(graph_filename + "_oslo_files", "tp"), G)
		return result

def clusterMOSES(G):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, "network.dat")
		output_filename = os.path.join(tempdir, "output.dat")
		graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListTabZero)
		subprocess.call([code_path + "/MOSES/moses-binary-linux-x86-64", graph_filename, output_filename])
		result = graphio.CoverReader().read(output_filename, G)
		return result

def clusterGCE(G, alpha = 1.5):
	with tempfile.TemporaryDirectory() as tempdir:
		graph_filename = os.path.join(tempdir, "network.edgelist")
		cover_filename = os.path.join(tempdir, "cover.txt")
		cover_file = open(cover_filename, 'x')
		graphio.writeGraph(G, graph_filename, fileformat=graphio.Format.EdgeListSpaceZero)
		subprocess.call([code_path + "/GCECommunityFinder/GCECommunityFinderUbuntu910", graph_filename, "4", "0.6", str(alpha), ".75"], stdout=cover_file)
		cover_file.close()
		C = graphio.CoverReader().read(cover_filename, G)
		# for u in G.nodes():
		# 	if not C.contains(u):
		# 		C.toSingleton(u)
		return C

def genLFR(N=1000, k=25, maxk=50, mu=0.01, t1=2, t2=1, minc=20, maxc=50, on=500, om=3, C=None):
	args = [code_path + "/lfr_graph_generator/benchmark", "-N", N, "-k", k, "-maxk", maxk, "-mu", mu, "-t1", t1, "-t2", t2, "-minc", minc, "-maxc", maxc]
	if on > 0:
		args.extend(["-on", on, "-om", om])
	if C is not None:
		args.extend(["-C", C])

	with tempfile.TemporaryDirectory() as tempdir:
		old_dir = os.getcwd()
		try:
			os.chdir(tempdir)
			with open("time_seed.dat", "w") as f:
				f.write(str(random.randint(0, 2**31)))
			subprocess.call(map(str, args), stdout=dev_null)
		finally:
			os.chdir(old_dir)

		G = graphio.readGraph(os.path.join(tempdir, "network.dat"), fileformat=graphio.Format.LFR)
		if on == 0:
			C = community.readCommunities(os.path.join(tempdir, "community.dat"), format='edgelist-s1')
		else:
			C = graphio.EdgeListCoverReader(1).read(os.path.join(tempdir, "community.dat"), G)
		return (G, C)

def getFacebookData(name, attribute):
	attribute_dict = {
		"student_fac" : 0,
		"gender" : 1,
		"major_index" : 2,
		"second_major" : 3,
		"dorm" : 4,
		"year" : 5,
		"high_school" : 6,
		}

	if attribute not in attribute_dict:
		raise Exception("Attribute {0} not found".format(attribute))

	fileName = home_path + "/graphs/facebook100/{0}.mat".format(name)
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

def getFacebookGraph(name):
	return graphio.readMat(home_path + "/graphs/facebook100/{0}.mat".format(name), key="A")

def getAmazonGraph():
	g = graphio.readGraph(code_path + "/graphs/com-amazon.ungraph.txt", fileformat=graphio.Format.SNAP)
	c = graphio.EdgeListCoverReader(1).read(code_path + "/graphs/com-amazon.all.dedup.cmty.txt", g)
	# c = community.readCommunities(code_path + "/graphs/com-amazon.all.dedup.cmty.txt", format='edgelist-t1')
	return g, c