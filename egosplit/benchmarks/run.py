from collections import OrderedDict

from networkit.community import PLM, PLP, LPPotts, OLP
from ..external import *
from .evaluation.metrics import write_results_to_file, add_compact_results, \
	print_compact_results
from .evaluation.ego_net_partition import analyse_ego_net_partitions
from .evaluation.stream_to_gephi import stream_partition
from .evaluation.cover_analysis import analyse_cover
from .cover_benchmark import CoverBenchmark
from .algorithms import *
from .graph import BenchGraph, LFRGraph


def start_benchmarks():
	iterations = 1
	append_results = False
	evaluations = [
		"metrics",
		"cover",
		# "ego_nets",
		# "stream_to_gephi",
	]
	store_ego_nets = "ego_nets" in evaluations \
	                 or "stream_to_gephi" in evaluations

	print("Creating Graphs...")
	graphs = get_graphs(iterations)
	algos = get_algos(store_ego_nets)

	print("Starting benchmarks...")
	result_summary = OrderedDict()
	if "stream_to_gephi" in evaluations:
		benchmarks = create_benchmarks(graphs, algos)
		run_benchmarks(benchmarks)
		evaluate_result(graphs, benchmarks, evaluations, append_results,
		                result_summary)
	else:
		for graph in graphs:
			for algo in algos:
				benchmarks = create_benchmarks([graph], [algo])
				run_benchmarks(benchmarks)
				evaluate_result(graphs, benchmarks, evaluations, append_results,
				                result_summary)
				append_results = True

	print_result_summary(result_summary)


def get_result_dir():
	return "./results/"


def print_result_summary(summary):
	if (summary):
		print_compact_results(summary, get_result_dir())


def evaluate_result(graphs, benchmarks, evaluations, append, summary):
	result_dir = get_result_dir()
	if "metrics" in evaluations:
		write_results_to_file(benchmarks, result_dir, append)
		compact_metrics = [
			"time",
			"nmi",
			"entropy",
			# "f1",
			# "f1_rev",
		]
		add_compact_results(summary, benchmarks, compact_metrics)
	if "cover" in evaluations:
		analyse_cover(graphs, benchmarks, result_dir, append)
	if "ego_nets" in evaluations:
		analyse_ego_net_partitions(benchmarks, result_dir, append)
	if "stream_to_gephi" in evaluations:
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
	for om in range(3, 6):
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
	algos.append(GroundTruth())
	# algos.append(OlpAlgorithm())
	# algos.append(GceAlgorithm())
	# algos.append(MosesAlgorithm())
	algos.append(OslomAlgorithm())

	partition_algos = OrderedDict()
	# partition_algos['PLP'] = [lambda g: PLP(g, 1, 20).run().getPartition()]
	# partition_algos['PLM_0.7'] = [lambda g: PLM(g, False, 0.7, "none").run().getPartition()]
	partition_algos['PLM_1.0'] = [lambda g: PLM(g, False, 1.0, "none").run().getPartition()]
	# partition_algos['PLM_1.2'] = [lambda g: PLM(g, False, 1.2, "none").run().getPartition()]
	# partition_algos['PLM_refine'] = [lambda g: PLM(g, True, 1.0, "none").run().getPartition()]
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
	algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="OSLOM")

	return algos


def get_ego_parameters(storeEgoNets):
	ego_parameters = OrderedDict()
	standard = {
		"weightFactor": 0,
		"weightOffset": 1,
		"storeEgoNet": "Yes" if storeEgoNets else "No",
		"addEgoNode": "No",
		"processEgoNet": "none",
		"addNodesFactor": 0,
		"addNodesExponent": 0,
	}
	extend_standard = {
		**standard,
		"addNodesFactor": 1,
		"addNodesExponent": 1,
		"processEgoNet": "extend",
		"edgesBetweenNeigNeig": "No",
		"extendRandom": "No",
	}
	edge_scores_standard = {
		**extend_standard,
		"extendStrategy": "edgeScore",
		"scoreStrategy": "score",
	}
	triangles_standard = {
		**extend_standard,
		"addNodesFactor": 10,
		"addNodesExponent": 0.6,
		"extendStrategy": "triangles",
		"scoreStrategy": "score_normed",
		"minTriangles": 2,
		"keepOnlyTriangles": "No",
		"edgesBetweenNeigNeig": "Yes",
		"triangleThreshold": 0,
	}

	ego_parameters['base'] = standard
	ego_parameters['edges'] = {
		**edge_scores_standard,
	}
	ego_parameters['triangles'] = {
		**triangles_standard,
	}

	# for score_strategy in ["score", "score_normed"]:
	# 	for add_nodes_factor in [2, 4, 6, 8, 10]:
	# 		ego_parameters['triangles_{}_{}'.format(add_nodes_factor, score_strategy)] = {
	# 			**triangles_standard,
	# 			"addNodesFactor": add_nodes_factor,
	# 			"scoreStrategy": score_strategy,
	# 		}
	# ego_parameters['triangles_noNN'] = {
	# 	**triangles_standard,
	# 	"edgesBetweenNeigNeig": "No",
	# }
	# add_nodes_factor_exponents = [
	# 	(10, 0.6),
	# 	(2, 1),
	# ]
	# for edges_NN in ["Yes"]:
	# 	for factor, exponent in add_nodes_factor_exponents:
	# 		name = 'triangles_{factor}*{exponent}{edges}'.format(
	# 			factor=factor,
	# 			exponent=exponent,
	# 			edges="_noNN" if edges_NN == "No" else ""
	# 		)
	# 		ego_parameters[name] = {
	# 			**triangles_standard,
	# 			"addNodesFactor": factor,
	# 			"addNodesExponent": exponent,
	# 			"edgesBetweenNeigNeig": edges_NN,
	# 		}

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

