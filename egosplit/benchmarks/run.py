from collections import OrderedDict

from networkit.community import PLM, PLP, LPPotts, OLP
from ..external import *
from .evaluation.metrics import write_results_to_file, print_benchmarks_compact
from .evaluation.ego_net_partition import analyse_ego_net_partitions
from .evaluation.stream_to_gephi import stream_partition
from .evaluation.cover_analysis import analyse_cover
from .cover_benchmark import CoverBenchmark
from .algorithms import *
from .graph import BenchGraph, LFRGraph


def start_benchmarks():
	iterations = 1
	store_ego_nets = False
	graphs = get_graphs(iterations)
	algos = get_algos(store_ego_nets)
	benchmarks = create_benchmarks(graphs, algos)

	run_benchmarks(benchmarks)

	result_dir = "./results/"
	append = False
	write_results_to_file(benchmarks, result_dir, append)
	compact_metrics = ["nmi"]
	print_benchmarks_compact(benchmarks, compact_metrics)
	analyse_cover(graphs, benchmarks, result_dir, append)
	if store_ego_nets:
		analyse_ego_net_partitions(benchmarks, result_dir, append)
		stream_partition(graphs, benchmarks)


def get_graphs(iterations):
	# ************************************************************************************
	# *                             Input graphs                                         *
	# ************************************************************************************
	graphs = []
	# graphs.append(BenchGraph(*getAmazonGraph5000(), "real_Amazon_5000"))
	# graphs.append(BenchGraph(*getAmazonGraph5000(True), "real_Amazon_5000_no_small"))
	# graphs.append(BenchGraph(*getAmazonGraphAll(), "real_Amazon_All"))
	# graphs.append(BenchGraph(*getAmazonGraphAll(True), "real_Amazon_All_no_small"))
	# graphs.append(BenchGraph(*getDBLPGraph(), "real_DBLP"))
	# graphs.append(BenchGraph(*getLiveJournalGraph(), 'real_LiveJournal'))
	# graphs.append(BenchGraph(*getOrkutGraph(), 'real_Orkut'))

	LFR_graph_args = OrderedDict()
	# LFR_graph_args['0.01'] = {'N': 1000, 'k': 25, 'maxk': 50, 'mu': 0.01, 'minc': 20,
	# 						'maxc': 50, 'on': 500, 'om': 3}
	# LFR_graph_args['0.1'] = {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.1, 'minc': 5,
	# 					   'maxc': 50, 'on': 100, 'om': 2}
	# LFR_graph_args['0.3'] = {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.3, 'minc': 5,
	# 					   'maxc': 50, 'on': 100, 'om': 2}
	for om in range(1, 6):
		LFR_graph_args['om_' + str(om)] = {
			'N': 2000, 'k': 18 * om, 'maxk': 120, 'minc': 60, 'maxc': 100,
			't1': 2, 't2': 2, 'mu': 0.2, 'on': 2000, 'om': om}
	# for mu_factor in range(0, 51, 5):
	# 	LFR_graph_args['mu_' + str(mu_factor).rjust(2,'0')] = {
	# 		'N': 1000, 'k': 20, 'maxk': 50, 'minc': 10, 'maxc': 50,
	# 		't1': 2, 't2': 1, 'mu': 0.01 * mu_factor, 'on': 1000, 'om': 2}
	create_LFR_graphs(graphs, LFR_graph_args, iterations)
	for graph in graphs:
		graph.graph.indexEdges()
	return graphs


def get_algos(storeEgoNets):
	# ************************************************************************************
	# *                         Benchmark algorithms                                     *
	# ************************************************************************************
	algos = []
	# algos.append(OslomAlgorithm())
	# algos.append(MosesAlgorithm())
	# algos.append(OlpAlgorithm())
	# algos.append(GceAlgorithm())

	partition_algos = OrderedDict()
	# partition_algos['PLP'] = [lambda g: PLP(g, 1, 20).run().getPartition()]
	partition_algos['PLM'] = [lambda g: PLM(g, False, 1.0, "none").run().getPartition()]
	# partition_algos['LPPotts'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition()]
	# partition_algos['LPPotts_par'] = [
	# 	lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition(),
	# 	lambda g: LPPotts(g, 0, 1, 20, True).run().getPartition()]
	# partition_algos['Infomap'] = [lambda g: clusterInfomap(g)]
	# partition_algos['Surprise'] = [lambda g: partitionLeiden(g, "surprise")]
	# partition_algos['Surprise_PLM'] = [
	# 	lambda g: partitionLeiden(g, "surprise"),
	# 	lambda g: PLM(g, False, 1.0, "none").run().getPartition()]

	ego_parameters = get_ego_parameters(storeEgoNets)

	algos += create_egosplit_algorithms(partition_algos, ego_parameters)
	# algos += create_egosplit_algorithms(partition_algos, clean_up="OSLOM")

	return algos


def get_ego_parameters(storeEgoNets):
	ego_parameters = OrderedDict()
	standard = {
		"weightFactor": "0",
		"weightOffset": "1",
		"addNodesFactor": "1",
		"addNodesExponent": "0.6",
		"storeEgoNet": "Yes" if storeEgoNets else "No",
		"addEgoNode": "No",
		"extendStrategy": "none",
	}
	ego_parameters['base'] = standard
	simple_nn_standard = {
		**standard,
		"extendStrategy": "simpleNN",
		"searchNeigOfNeighInDirectedGraph": "No",
		"discardNeigOfNeigEdgesAtFirst": "Yes",
		"discardNonTriangle": "Yes",
		"removeLowDegreeNeig": "Yes",
		"minDegreeCleaning": "1",
		"edgesBetweenNeigNeig": "No",
	}
	edge_scores_standard = {
		**standard,
		"extendStrategy": "edgeScores",
		"discardThreshold": "0.0",
		"scoreStrategy": "count",
	}
	ego_parameters['base_0.1'] = {
		**standard,
		"weightFactor": "1",
		"weightOffset": "0",
		"discardThreshold": 0.1,
	}
	ego_parameters['base_0.2'] = {
		**standard,
		"weightFactor": "1",
		"weightOffset": "0",
		"discardThreshold": 0.2,
	}
	ego_parameters["edgeScores"] = edge_scores_standard
	ego_parameters["edgeScores_0.1"] = {
		**edge_scores_standard,
		"weightFactor": "1",
		"weightOffset": "0",
		"discardThreshold": 0.1,
	}
	ego_parameters["edgeScores_0.2"] = {
		**edge_scores_standard,
		"weightFactor": "1",
		"weightOffset": "0",
		"discardThreshold": 0.2,
	}
	# ego_parameters["edgeScores_weighted_geometric"] = {
	# 	**edge_scores_standard,
	# 	"weightFactor": "1",
	# 	"weightOffset": "0",
	# 	"scoreStrategy": "geometric",
	# }
	# ego_parameters["edgeScores_weighted_numScores"] = {
	# 	**edge_scores_standard,
	# 	"weightFactor": "1",
	# 	"weightOffset": "0",
	# 	"scoreStrategy": "numScores",
	# }
	# add_nodes_exponent_factors = {
	# 	"0.4": 8,
	# 	"0.6": 5,
	# 	"0.8": 3,
	# 	"1.0": 2,
	# }
	# for exponent, max_factor in add_nodes_exponent_factors.items():
	# 	num_steps = 20
	# 	step_size = max_factor / num_steps
	# 	for i in range(0, num_steps):
	# 		factor = str(i * step_size)[:5]
	# 		ego_parameters['edgeScore_' + exponent + "_" + factor] = {
	# 			**edge_scores_standard,
	# 			"addNodesExponent": exponent,
	# 			"addNodesFactor": factor,
	# 		}

	# for i in range(4, 11, 2):
	# 	value = str(0.1 * i)[:5]

	return ego_parameters


def create_LFR_graphs(graphs, graphArgs, iterations):
	for argsName in graphArgs.keys():
		args = graphArgs[argsName]
		graphName = 'LFR_' + argsName
		for i in range(iterations):
			graph_wrapper = LFRGraph(graphName, args)
			graphs.append(graph_wrapper)
	return graphs


def create_egosplit_algorithms(partition_algos, ego_parameters, clean_up=""):
	algos = []
	for part_name in partition_algos:
		for para_name, parameters in ego_parameters.items():
			algos.append(EgoSplitAlgorithm(
				part_name + "_" + para_name,
				{key.encode("utf-8"): str(value).encode("utf-8")
				 for (key, value) in parameters.items()},
				*partition_algos[part_name],
				clean_up=clean_up))
	return algos


def create_benchmarks(graphs, algos):
	benchmarks = []
	for algo in algos:
		for graph in graphs:
			benchmarks.append(CoverBenchmark(algo.copy(), graph))
	return benchmarks


def run_benchmarks(benchmarks):
	for benchmark in benchmarks:
		benchmark.run()

