from egosplit.external import *
from networkit.community import EgoSplitting, CoverF1Similarity, PLM, PLP, \
	LPPotts, OLP
from networkit import setLogLevel, setNumberOfThreads
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


def print_results_to_file(results, file, print_head):
	algos = list(results.keys())
	graphs = list(results[algos[0]].keys())
	metrics = results[algos[0]][graphs[0]].keys()
	result_decimals = 4
	str_width = max(4 + result_decimals, *[len(metric) for metric in metrics]) + 1
	algo_width = max([len(algo) for algo in algos]) + 2
	graph_width = max([len(graph) for graph in graphs]) + 2

	if print_head:
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
	result_decimals = 4
	str_width = max(4 + result_decimals, *[len(metric) for metric in metrics]) + 3
	str_width_first = max([len(algo) for algo in algos]) + 1

	# Header
	graph_header = "\n" + "Graph ".ljust(str_width_first)
	for graph in graphs:
		graph_header += "| " + graph.rjust(len(metrics) * str_width - 3) + ' '
	graph_header += '\n' + str().ljust(str_width_first
									   + len(metrics) * str_width * len(graphs), '-')
	graph_header += '\n' + str().ljust(str_width_first)
	for _ in graphs:
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


def run_benchmarks(algos, graphs):
	results = OrderedDict()
	for algo in algos:
		name = algo.getName()
		print("\t" + name)
		results[name] = OrderedDict()
		for graphName in graphs:
			print("\t\t" + graphName)
			results[name][graphName] = OrderedDict()
			f1Sum = 0
			f1_rev_sum = 0
			nmi_sum = 0
			timeSum = 0
			iterations = len(graphs[graphName])
			for graph, cover in graphs[graphName]:
				algo.run(graph)
				f1 = calc_F1(graph, algo.getCover(), cover)
				f1_rev = calc_F1(graph, cover, algo.getCover())
				nmi = calc_NMI(graph, algo.getCover(), cover)
				t = algo.getTime()
				f1Sum += f1
				f1_rev_sum += f1_rev
				nmi_sum += nmi
				timeSum += t
				print("\t\t\tF1:" + str(f1) + ", NMI:" + str(nmi) + ", t:" + str(t))
			results[name][graphName]['F1'] = f1Sum / iterations
			results[name][graphName]['F1_(rev)'] = f1_rev_sum / iterations
			results[name][graphName]['NMI'] = nmi_sum / iterations
			results[name][graphName]['time'] = timeSum / iterations
	return results


def start_benchmarks():
	# setNumberOfThreads(1)
	iterations = 20
	graphs = OrderedDict()
	LFR_graph_args = OrderedDict()
	algos = []
	partition_algos = OrderedDict()

	append = True
	open_mode = 'w'
	if append:
		open_mode = 'a'
	ego_file = open('ego_timings.txt', open_mode)
	if not append:
		s = 21
		ego_file.write('create EgoNets'.ljust(s) + 'a) assign IDs'.ljust(s)
					   + 'b) find triangles'.ljust(s) + 'c) cluster EgoNet'.ljust(s)
					   + 'd) compact EgoNet'.ljust(s) + 'e) EgoNet subsets'.ljust(s)
					   + 'f) EgoNet subsetCnt'.ljust(s) + 'g) clean up'.ljust(s)
					   + 'split Personas'.ljust(s) + 'connect Personas'.ljust(s)
					   + 'Persona Clustering'.ljust(s) + 'create Cover'.ljust(s)
					   + '\n')

	# ************************************************************************************
	# * Input graphs                                                                     *
	# ************************************************************************************
	# graphs['Amazon'] = [getAmazonGraph()]
	# graphs['DBLP'] = [getDBLPGraph()]
	# graphs['LiveJournal'] = [getLiveJournalGraph()]
	# graphs['Orkut'] = [getOrkutGraph()]
	LFR_graph_args['0.01'] = {'N': 1000, 'k': 25, 'maxk': 50, 'mu': 0.01, 'minc': 20,
							'maxc': 50, 'on': 500, 'om': 3}
	LFR_graph_args['0.1'] = {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.1, 'minc': 5,
						   'maxc': 50, 'on': 100, 'om': 2}
	LFR_graph_args['0.3'] = {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.3, 'minc': 5,
						   'maxc': 50, 'on': 100, 'om': 2}
	create_LFR_graphs(graphs, LFR_graph_args, iterations)

	# ************************************************************************************
	# * Benchmark algorithms                                                             *
	# ************************************************************************************
	partition_algos['PLP'] = [lambda g: PLP(g, 1, 20).run().getPartition()]
	partition_algos['PLM'] = [lambda g: PLM(g, False, 1.0, "none").run().getPartition()]
	partition_algos['LPPotts'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition()]
	partition_algos['LPPotts_par'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition(),
									  lambda g: LPPotts(g, 0, 1, 20, True).run().getPartition()]
	for partitionAlgo in partition_algos:
		algos.append(EgoSplitAlgorithm(ego_file, partitionAlgo, *partition_algos[partitionAlgo]))
	algos.append(OLPAlgorithm())
	algos.append(GCEAlgorithm())
	# algos.append(MOSESAlgorithm())
	# algos.append(OSLOMAlgorithm())

	# ************************************************************************************
	# * Run benchmarks and print results                                                 *
	# ************************************************************************************
	results = run_benchmarks(algos, graphs)

	result_file = open('results.txt', open_mode)
	print_results_to_file(results, result_file, not append)
	print_results_compact(results)

# setLogLevel("INFO")
start_benchmarks()