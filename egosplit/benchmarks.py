from egosplit.external import *
from networkit.community import EgoSplitting, CoverF1Similarity, PLM, PLP, \
	LPPotts, OLP
from networkit import setLogLevel
import collections


def calcSimilarity(algo, graph, refCover):
	cover = algo(graph)
	similarity = CoverF1Similarity(graph, cover, refCover).run()
	return similarity.getUnweightedAverage()


def createLFRGraphs(graphs, graphArgs, iterations):
	for argsName in graphArgs.keys():
		args = graphArgs[argsName]
		graphName = 'LFR_' + argsName
		graphs[graphName] = []
		for i in range(0, iterations):
			graphs[graphName].append(genLFR(**args))
	return graphs


def createEgoAlgos(partitionAlgos):
	algos = []

	def buildEgoLambda(partAlgoName):
		return lambda g: EgoSplitting(g, partitionAlgos[partAlgoName][0],
									  partitionAlgos[partAlgoName][1]) \
			.run().getCover()

	for partAlgoName in sorted(partitionAlgos.keys()):
		algos.append(("EgoSplit(" + partAlgoName + ")",
					  buildEgoLambda(partAlgoName)))
	return algos


def runBenchmarks(graphs, algos, iterations):
	results = collections.OrderedDict()
	for graphName in graphs.keys():
		results[graphName] = collections.OrderedDict()
		print(graphName)
		for algo in algos:
			print("\t" + algo[0])
			simSum = 0
			for graph, cover in graphs[graphName]:
				s = calcSimilarity(algo[1], graph, cover)
				simSum += s
				print("\t\t" + str(s))
			avg = simSum / iterations
			results[graphName][algo[0]] = avg
	return results


def printResults(results, algos):
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


def startBenchmarks():
	# Input graphs
	graphs = collections.OrderedDict()
	graphs['Amazon'] = [getAmazonGraph()]
	# graphs['DBLP'] = [getDBLPGraph()]
	# graphs['LiveJournal'] = [getLiveJournalGraph()]
	# graphs['Orkut'] = [getOrkutGraph()]
	LFRgraphArgs = collections.OrderedDict()
	# LFRgraphArgs['0.01'] = {'N': 1000, 'k': 25, 'maxk': 50, 'mu': 0.01, 'minc': 20,
	# 						'maxc': 50, 'on': 500, 'om': 3}
	# LFRgraphArgs['0.1'] = {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.1, 'minc': 5,
	# 					   'maxc': 50, 'on': 100, 'om': 2}
	# LFRgraphArgs['0.3'] = {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.3, 'minc': 5,
	# 					   'maxc': 50, 'on': 100, 'om': 2}

	# Benchmark algorithms
	partitionAlgos = collections.OrderedDict({
		# 'PLP': [lambda g: PLP(g).run().getPartition(), None],
		# 'PLM': [lambda g: PLM(g).run().getPartition(), None],
		# 'LPPotts': [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition(), None],
		# 'LPPotts_par': [lambda g: LPPotts(g, 0.1, 1, 50, True).run().getPartition(), None],
	})
	algos = [
		# ("OSLOM", lambda g: clusterOSLOM(g)),
		# ("GCE", lambda g: clusterGCE(g)),
		("OLP", lambda g: OLP(g).run().getCover()),
		("MOSES", lambda g: clusterMOSES(g)),
	]

	iterations = 1

	# Create input graphs and benchmark algorithms
	createLFRGraphs(graphs, LFRgraphArgs, iterations)
	algos.extend(createEgoAlgos(partitionAlgos))

	results = runBenchmarks(graphs, algos, iterations)
	printResults(results, algos)

startBenchmarks()
