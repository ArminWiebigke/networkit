from egosplit.external import genLFR, clusterOSLOM, clusterGCE, clusterMOSES
from networkit.community import EgoSplitting, CoverF1Similarity, PLM, PLP, \
	LPPotts
import collections


def calcSimilarity(algo, graph, refCover):
	cover = algo(graph)
	similarity = CoverF1Similarity(graph, cover, refCover).run()
	return similarity.getUnweightedAverage()


def createGraphs(graphNames, graphArgs, iterations):
	graphs = {}
	for graphName in graphNames:
		args = graphArgs[graphName]
		graphs[graphName] = []
		for i in range(0, iterations):
			graphs[graphName].append(genLFR(**args))
	return graphs


def createAlgos(partitionAlgos, partitionAlgoNames):
	algos = [
		# ("OSLOM", lambda g: clusterOSLOM(g)),
		# ("GCE", lambda g: clusterGCE(g)),
		# ("MOSES", lambda g: clusterMOSES(g))
	]

	def buildEgoLambda(partAlgoName):
		return lambda g: EgoSplitting(g, partitionAlgos[partAlgoName][0],
									  partitionAlgos[partAlgoName][1]) \
			.run().getCover()

	for partAlgoName in partitionAlgoNames:
		algos.append(("EgoSplit(" + partAlgoName + ")",
					  buildEgoLambda(partAlgoName)))
	return algos


def startBenchmarks():
	# Benchmark parameters
	graphArgs = {
		'0.01': {'N': 1000, 'k': 25, 'maxk': 50, 'mu': 0.01, 'minc': 20,
				 'maxc': 50, 'on': 500, 'om': 3},
		'0.1': {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.1, 'minc': 5,
				'maxc': 50, 'on': 100, 'om': 2},
		'0.3': {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.3, 'minc': 5,
				'maxc': 50, 'on': 100, 'om': 2},
	}
	partitionAlgos = collections.OrderedDict({
		'PLP': [lambda g: PLP(g).run().getPartition(), None],
		'PLM': [lambda g: PLM(g).run().getPartition(), None],
		'LPPotts': [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition(), None],
	})

	for a in range(1, 6):
		for k in [5, 10, 20]:
			partitionAlgos['LPPotts_' + str(a) + "_" + str(k)] = [lambda g, a=a, k=k: LPPotts(g, 0.1, a, k).run().getPartition(), None]


	graphNames = [
		'0.01',
		'0.1',
		'0.3'
	]
	partitionAlgoNames = [
		'PLP',
		'PLM',
		'LPPotts',
	]
	iterations = 30

	# Create input graphs and benchmark algorithms
	graphs = createGraphs(graphNames, graphArgs, iterations)
	algos = createAlgos(partitionAlgos, partitionAlgos.keys())

	# Run benchmarks
	results = collections.OrderedDict()
	for graphName in graphNames:
		results[graphName] = collections.OrderedDict()
		for algo in algos:
			print("\t" + algo[0])
			simSum = 0
			for graph, cover in graphs[graphName]:
				s = calcSimilarity(algo[1], graph, cover)
				simSum += s
				print("\t\t" + str(s))
			avg = simSum / iterations
			results[graphName][algo[0]] = avg

	# Print results
	str_width = 22
	str_width_first = 25
	graph_header = "\n" + "Graph ".ljust(str_width_first)
	for graph in results.keys():
		graph_header += ("| " + graph).ljust(str_width)
	print(graph_header)
	print(str().ljust(str_width_first + str_width * len(results.keys()), '-'))
	for algo, _ in algos:
		algo_results = (algo + " ").ljust(str_width_first)
		for graph in results.keys():
			algo_results += ("| " + str(results[graph][algo])).ljust(str_width)
		print(algo_results)


startBenchmarks()
