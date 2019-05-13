from collections import OrderedDict

from networkit.community import PLM, PLP, LPPotts, OLP
from ..external import *
from .evaluation.metrics import write_results_to_file, add_compact_results, \
	print_compact_results
from .evaluation.ego_net_partition import analyse_ego_net_partitions
from .evaluation.stream_to_gephi import stream_partition
from .evaluation.cover_analysis import analyse_cover
from .evaluation.cleanup import cleanup_test
from .cover_benchmark import CoverBenchmark, MetricCache
from .algorithms import *
from .graph import BenchGraph, LFRGraph
import egosplit.benchmarks.evaluation.benchmark_metric as bm


def start_benchmarks():
	iterations = 1
	append_results = False
	#append_results = True
	evaluations = [
		"metrics",
		"cover",
		"ego_nets",
		# "stream_to_gephi",
		# "cleanup",
	]
	stream_to_gephi = "stream_to_gephi" in evaluations
	store_ego_nets = "ego_nets" in evaluations \
	                 or stream_to_gephi

	print("Creating Graphs...")
	graphs = get_graphs(iterations)
	algos = get_algos(store_ego_nets)
	if stream_to_gephi and len(graphs) * len(algos) > 8:
		raise RuntimeError("Too many runs to stream!")

	print("Starting benchmarks...")
	result_summary = OrderedDict()
	if stream_to_gephi:
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
	if summary:
		print_compact_results(summary, get_result_dir())


# TODO: Mehr output: GraphID, timestamp
def evaluate_result(graphs, benchmarks, evaluations, append, summary):
	result_dir = get_result_dir()
	if "metrics" in evaluations:
		metric_caches = [MetricCache(b) for b in benchmarks]
		metrics = [
			bm.Time,
			bm.NMI,
			# bm.Entropy,
			# bm.Entropy2,
			# bm.Entropy3,
			# bm.Entropy4,
			bm.F1,
			bm.F1_rev,
		]
		write_results_to_file(metric_caches, result_dir, metrics, append)
		compact_metrics = [
			bm.Time,
			bm.NMI,
			# bm.Entropy,
			# bm.Entropy2,
			# bm.Entropy3,
			# bm.Entropy4,
			# bm.F1,
			# bm.F1_rev,
		]
		add_compact_results(summary, metric_caches, compact_metrics)
	if "cover" in evaluations:
		analyse_cover(benchmarks, result_dir, append)
	if "ego_nets" in evaluations:
		analyse_ego_net_partitions(benchmarks, result_dir, append)
	if "stream_to_gephi" in evaluations:
		stream_partition(graphs, benchmarks)
	if "cleanup" in evaluations:
		cleanup_test(benchmarks, result_dir)


# ************************************************************************************
# *                             Input graphs                                         *
# ************************************************************************************
def get_graphs(iterations):
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
	for om in range(5, 6):
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


# ************************************************************************************
# *                         Benchmark algorithms                                     *
# ************************************************************************************
def get_algos(storeEgoNets):
	algos = []
	algos.append(GroundTruth())
	# algos.append(OlpAlgorithm())
	# algos.append(GceAlgorithm())
	# algos.append(MosesAlgorithm())
	# algos.append(OslomAlgorithm())

	partition_algos = OrderedDict()
	# partition_algos['PLP'] = [lambda g: PLP(g, 1, 20).run().getPartition()]
	# partition_algos['PLM_1.0'] = [lambda g: PLM(g, True, 1.0, "none").run().getPartition()]
	partition_algos['PLM_Info'] = [
		lambda g: PLM(g, True, 1.0, "none").run().getPartition(),
		lambda g: clusterInfomap(g)]
	# partition_algos['PLM_refine'] = [lambda g: PLM(g, True, 1.0, "none").run().getPartition()]
	# partition_algos['Potts'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition()]
	# partition_algos['LPPotts_par'] = [
	# 	lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition(),
	# 	lambda g: LPPotts(g, 0, 1, 20, True).run().getPartition()]
	# partition_algos['Infomap'] = [lambda g: clusterInfomap(g)]
	# partition_algos['Surprise'] = [lambda g: partitionLeiden(g, "surprise")]
	# partition_algos['Leiden_Mod'] = [lambda g: partitionLeiden(g, "modularity")]

	# for p_algos in partition_algos.values():
	# 	p_algos.append(lambda g: PLM(g, True, 1.0, "none").run().getPartition())
	ego_parameters = get_ego_parameters(storeEgoNets)

	algos += create_egosplit_algorithms(partition_algos, ego_parameters)
	algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="OSLOM")
	algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="OSLOM_large")
	algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="OSLOM_keep_1")
	algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="OSLOM_keep_10")
	algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="OSLOM_keep_25")
	algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="OSLOM_merge")
	algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="trim")
	#algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="trim_overl")
	algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="trim_merge_gt")
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="trim_merge")
	#algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="trim_merge_all")
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="trim_merge_overl")
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="tr_ol_mrg_ol")
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="trim_remove")
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="trim_ol_remove")
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="entropy")
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="entropy_merge")
	#algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up="entropy_all")

	return algos


def get_ego_parameters(store_ego_nets):
	ego_parameters = OrderedDict()
	standard = {
		"weightFactor": 0,
		"weightOffset": 1,
		"storeEgoNet": "Yes" if store_ego_nets else "No",
		"addEgoNode": "No",
		"processEgoNet": "none",
		"addNodesFactor": 0,
		"addNodesExponent": 0,
		"weightedEgoNet": "No",
		"partitionFromGroundTruth": "No",
	}
	extend_standard = {
		**standard,
		"addNodesFactor": 4,
		"addNodesExponent": 0.6,
		"processEgoNet": "extend",
		"edgesBetweenNeigNeig": "Yes",
		"extendRandom": "No",
		"minNodeDegree": 0,
		"triangleThreshold": 0,
		"extendOverDirected": "No",
		"keepOnlyTriangles": "No",
	}
	edge_scores_standard = {
		**extend_standard,
		"extendStrategy": "edgeScore",
		"scoreStrategy": "score^2_normed",
	}
	triangles_standard = {
		**extend_standard,
		"extendStrategy": "triangles",
		"scoreStrategy": "score_normed",
		"edgesBetweenNeigNeig": "Yes",
		"minNodeDegree": 0,
		"minTriangles": 0,
		"triangleThreshold": 0,
	}

	# ego_parameters['base'] = standard
	# ego_parameters['gt'] = {
	# 	**standard,
	# 	"partitionFromGroundTruth": "Yes",
	# }
	ego_parameters[''] = {
		**edge_scores_standard,
	}
	# ego_parameters['edges_dir'] = {
	# 	**edge_scores_standard,
	# 	"extendOverDirected": "Yes",
	# }
	# ego_parameters['triangles'] = {
	# 	**triangles_standard,
	# }
	# for threshold in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1, 1.2, 1.4, 1.7, 2]:
	# 	ego_parameters['triangles_{}'.format(threshold)] = {
	# 		**triangles_standard,
	# 		"triangleThreshold": threshold,
	# 	}

	# add_nodes_factor_exponents = [
	# 	(0, 0.6),
	# 	(1, 0.6),
	# 	(2, 0.6),
	# 	(3, 0.6),
	# 	(4, 0.6),
	# 	(6, 0.6),
	# 	(8, 0.6),
	# 	(12, 0.6),
	# 	(16, 0.6),
	# 	(0, 0.8),
	# 	(0.5, 0.8),
	# 	(1, 0.8),
	# 	(1.5, 0.8),
	# 	(2, 0.8),
	# 	(3, 0.8),
	# 	(4, 0.8),
	# 	(6, 0.8),
	# 	(8, 0.8),
	# 	(0, 1),
	# 	(0.25, 1),
	# 	(0.5, 1),
	# 	(0.75, 1),
	# 	(1, 1),
	# 	(1.5, 1),
	# 	(2, 1),
	# 	(2.5, 1),
	# 	(3, 1),
	# 	(4, 1),
	# ]
	# for factor, exponent in add_nodes_factor_exponents:
	# 	name = 'edges_f-{factor}*e-{exponent}'.format(
	# 		factor=factor,
	# 		exponent=exponent,
	# 	)
	# 	ego_parameters[name] = {
	# 		**edge_scores_standard,
	# 		"addNodesFactor": factor,
	# 		"addNodesExponent": exponent,
	# 	}

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
			name = "{}{}{}".format(part_name, "_" if para_name else "", para_name)
			algos.append(EgoSplitAlgorithm(
				name,
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
