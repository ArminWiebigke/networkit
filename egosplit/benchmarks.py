from egosplit.external import *
from networkit.community import EgoSplitting, CoverF1Similarity, PLM, PLP, \
	LPPotts, OLP
from networkit import setLogLevel
from collections import OrderedDict
from egosplit.algorithms import *

def calc_F1(graph, cover, refCover):
	similarity = CoverF1Similarity(graph, cover, refCover).run()
	return similarity.getUnweightedAverage()


def create_LFR_graphs(graphs, graphArgs, iterations):
	for argsName in graphArgs.keys():
		args = graphArgs[argsName]
		graphName = 'LFR_' + argsName
		graphs[graphName] = []
		for i in range(0, iterations):
			graphs[graphName].append(genLFR(**args))
	return graphs


def fixed_decimals(val, decimals = 3):
	string = str(val)
	split_strs = string.split('.')
	split_strs[1] = split_strs[1].ljust(decimals, '0')
	return split_strs[0] + '.' + split_strs[1][:decimals]


def print_results_to_file(results):
	algos = list(results.keys())
	graphs = list(results[algos[0]].keys())
	metrics = results[algos[0]][graphs[0]].keys()
	result_decimals = 3
	str_width = 6 + result_decimals
	algo_width = max([len(algo) for algo in algos]) + 2
	graph_width = max([len(graph) for graph in graphs]) + 2

	file = open('results.txt', 'w')
	header = 'algo'.ljust(algo_width) + 'graph'.ljust(graph_width)
	for metric in metrics:
		header += metric.ljust(str_width)
	file.write(header + '\n')
	print(header)
	for algo in algos:
		for graph in graphs:
			result_line = algo.ljust(algo_width) + graph.ljust(graph_width)
			for metric in metrics:
				result_line += fixed_decimals(results[algo][graph][metric],
											  result_decimals).ljust(str_width)
			file.write(result_line + '\n')
			print(result_line)


def print_results_compact(results):
	algos = list(results.keys())
	graphs = list(results[algos[0]].keys())
	metrics = results[algos[0]][graphs[0]].keys()
	result_decimals = 3
	str_width = 7 + result_decimals
	str_width_first = max([len(algo) for algo in algos]) + 1

	# Header
	graph_header = "\n" + "Graph ".ljust(str_width_first)
	for graph in graphs:
		graph_header += "| " + graph.rjust(len(metrics) * str_width - 3) + ' '
	graph_header += '\n' + str().ljust(str_width_first
									   + len(metrics) * str_width * len(graphs), '-')
	graph_header += '\n' + str().ljust(str_width_first)
	for graph in graphs:
		for metric in metrics:
			graph_header += "| " + metric.rjust(str_width - 3) + ' '
	graph_header += '\n' + str().ljust(
		str_width_first + len(metrics) * str_width * len(graphs), '-')
	print(graph_header)

	# Results
	for algo in algos:
		algo_results = (algo + " ").ljust(str_width_first)
		for graph in graphs:
			for metric in metrics:
				algo_results += "| " + str(
					fixed_decimals(results[algo][graph][metric],
								   result_decimals)).rjust(str_width - 3) + ' '
		print(algo_results)


def run_benchmarks(algos, graphs, iterations):
	results = OrderedDict()
	for algo in algos:
		name = algo.getName()
		print("\t" + name)
		results[name] = OrderedDict()
		for graphName in graphs:
			print("\t\t" + graphName)
			results[name][graphName] = OrderedDict()
			f1Sum = 0
			nmiSum = 0
			nmiSum2 = 0
			timeSum = 0
			for graph, cover in graphs[graphName]:
				algo.run(graph)
				f1 = calc_F1(graph, algo.getCover(), cover)
				nmi, nmi2 = calcNMI(graph, algo.getCover(), cover)
				t = algo.getTime()
				f1Sum += f1
				nmiSum += nmi
				nmiSum2 += nmi2
				timeSum += t
				print("\t\t\tF1:" + str(f1) + ", NMI:" + str(nmi) + ", t:" + str(t))
			results[name][graphName]['F1'] = f1Sum / iterations
			results[name][graphName]['NMI'] = nmiSum / iterations
			results[name][graphName]['NMI(2)'] = nmiSum2 / iterations
			results[name][graphName]['time'] = timeSum / iterations
	return results


def start_benchmarks():
	iterations = 1
	graphs = OrderedDict()
	LFRgraphArgs = OrderedDict()
	algos = []
	partitionAlgos = OrderedDict()

	# ************************************************************************************
	# * Input graphs                                                                     *
	# ************************************************************************************
	graphs['Amazon'] = [getAmazonGraph()]
	# graphs['DBLP'] = [getDBLPGraph()]
	# graphs['LiveJournal'] = [getLiveJournalGraph()]
	# graphs['Orkut'] = [getOrkutGraph()]
	# LFRgraphArgs['0.01'] = {'N': 1000, 'k': 25, 'maxk': 50, 'mu': 0.01, 'minc': 20,
	# 						'maxc': 50, 'on': 500, 'om': 3}
	# LFRgraphArgs['0.1'] = {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.1, 'minc': 5,
	# 					   'maxc': 50, 'on': 100, 'om': 2}
	# LFRgraphArgs['0.3'] = {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.3, 'minc': 5,
	# 					   'maxc': 50, 'on': 100, 'om': 2}
	create_LFR_graphs(graphs, LFRgraphArgs, iterations)

	# ************************************************************************************
	# * Benchmark algorithms                                                             *
	# ************************************************************************************
	partitionAlgos['PLP'] = [lambda g: PLP(g, 1, 20).run().getPartition()]
	partitionAlgos['PLM'] = [lambda g: PLM(g).run().getPartition()]
	# partitionAlgos['LPPotts'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition()]
	# partitionAlgos['LPPotts_par'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition(),
	# 								 lambda g: LPPotts(g, 0, 1, 20, True).run().getPartition()]
	for partitionAlgo in partitionAlgos:
		algos.append(EgoSplitAlgorithm(partitionAlgo, *partitionAlgos[partitionAlgo]))
	# algos.append(OLPAlgorithm())
	# algos.append(MOSESAlgorithm())
	# algos.append(GCEAlgorithm())
	# algos.append(OSLOMAlgorithm())

	# ************************************************************************************
	# * Run benchmarks and print results                                                 *
	# ************************************************************************************
	results = run_benchmarks(algos, graphs, iterations)
	print_results_to_file(results)
	print_results_compact(results)

# setLogLevel("INFO")
start_benchmarks()