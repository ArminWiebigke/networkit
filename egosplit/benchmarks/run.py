from collections import OrderedDict
from copy import copy

import egosplit.benchmarks.evaluation.benchmark_metric as bm

from networkit.stopwatch import clockit
from networkit.community import PLM, PLP, LPPotts, SLPA, OslomCleanUp
from networkit.structures import Cover
from networkit import setLogLevel
from ..external import *
from .evaluation.metrics import write_results_to_file, add_compact_results, \
	print_compact_results
from .evaluation.ego_net_partition import analyse_ego_net_partitions
from .evaluation.stream_to_gephi import stream_partition
from .evaluation.cover_analysis import analyse_cover
from egosplit.benchmarks.cleanup import cleanup_test
from .cover_benchmark import CoverBenchmark, MetricCache
from .algorithms import *
from .graph import BenchGraph, LFRGraph


def start_benchmarks():
	# setLogLevel('INFO')
	iterations = 5
	append_results = False
	# append_results = True
	evaluations = [
		'metrics',
		'cover',
		'ego_nets',
		# 'stream_to_gephi',
		# 'cleanup',
	]
	stream_to_gephi = 'stream_to_gephi' in evaluations
	store_ego_nets = 'ego_nets' in evaluations \
	                 or stream_to_gephi

	print('Creating Graphs...')
	graphs = get_graphs(iterations)
	algos = get_algos(store_ego_nets)
	if stream_to_gephi and len(graphs) * len(algos) > 8:
		raise RuntimeError('Too many runs to stream!')

	print('Starting benchmarks...')
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
	return './results/'


def print_result_summary(summary):
	if summary:
		print_compact_results(summary, get_result_dir())


# TODO: Mehr output: GraphID, timestamp
@clockit
def evaluate_result(graphs, benchmarks, evaluations, append, summary):
	result_dir = get_result_dir()
	if 'metrics' in evaluations:
		metric_caches = [MetricCache(b) for b in benchmarks]
		metrics = [
			bm.Time,
			bm.NMI,
			bm.F1,
			bm.F1_rev,
			# bm.Entropy,
			# bm.Entropy2,
			# bm.Entropy3,
			# bm.Entropy4,
		]
		write_results_to_file(metric_caches, result_dir, metrics, append)
		compact_metrics = [
			bm.Time,
			bm.NMI,
			bm.F1,
			bm.F1_rev,
			# bm.Entropy,
			# bm.Entropy2,
			# bm.Entropy3,
			# bm.Entropy4,
		]
		add_compact_results(summary, metric_caches, compact_metrics)
	if 'cover' in evaluations:
		analyse_cover(benchmarks, result_dir, append)
	if 'ego_nets' in evaluations:
		analyse_ego_net_partitions(benchmarks, result_dir, append)
	if 'stream_to_gephi' in evaluations:
		stream_partition(graphs, benchmarks)
	if 'cleanup' in evaluations:
		cleanup_test(benchmarks, result_dir)


# ****************************************************************************************
# *                                  Input graphs                                        *
# ****************************************************************************************
def get_graphs(iterations):
	graphs = []
	# graphs.append(BenchGraph(*getAmazonGraph5000(), 'Amazon_5000'))
	# graphs.append(BenchGraph(*getAmazonGraph5000(True), 'Amazon_5000_no_small'))
	# graphs.append(BenchGraph(*getAmazonGraphAll(), 'Amazon_All'))
	# graphs.append(BenchGraph(*getAmazonGraphAll(True), 'Amazon_All_no_small'))
	# graphs.append(BenchGraph(*getDBLPGraph(), 'DBLP'))
	# graphs.append(BenchGraph(*getLiveJournalGraph(), 'LiveJournal'))
	# graphs.append(BenchGraph(*getOrkutGraph(), 'Orkut'))
	# graphs.append(BenchGraph(*getFacebookGraph('Caltech36', clean=True), 'Caltech36'))
	# graphs.append(BenchGraph(*getFacebookGraph('Rice31', clean=True), 'Rice31'))
	# graphs.append(BenchGraph(*getFacebookGraph('Auburn71', clean=True), 'Auburn71'))
	# graphs.append(BenchGraph(*getFacebookGraph('Smith60', clean=True), 'Smith60'))
	# graphs.append(BenchGraph(*getFacebookGraph('Oklahoma97', clean=True), 'Oklahoma97'))

	LFR_graph_args = OrderedDict()

	# TODO: Graphs with different community sizes, e.g. 10-100
	for minc in [60]:
		for om in range(5, 6):
			for mu in [0.2]:
				k = 15 * om  # Number of neighbors per community independent of mixing factor
				k /= (1 - mu)
				maxk = 1.5 * k  # Scale max degree with average degree
				LFR_graph_args['mu_{}_om_{}'.format(mu, om)] = {
					'N': 2000, 'k': k, 'maxk': maxk, 'minc': minc, 'maxc': 100,
					't1': 2, 't2': 2, 'mu': mu, 'on': 2000, 'om': om}

	# # GCE paper
	# for om in range(1, 6):
	# 	LFR_graph_args['om_' + str(om)] = {
	# 		'N': 2000, 'k': 18 * om, 'maxk': 120, 'minc': 60, 'maxc': 100,
	# 		't1': 2, 't2': 2, 'mu': 0.2, 'on': 2000, 'om': om}
	# for om in range(1, 7):
	# 	k = 12 * om
	# 	LFR_graph_args['om_{}'.format(om)] = {
	# 		'N': 2000, 'k': k, 'maxk': 1.5 * k, 'minc': 60, 'maxc': 100,
	# 		't1': 2, 't2': 2, 'mu': 0.2, 'on': 1000, 'om': om}
	#
	# # Overlapping CD study
	# small = (10, 50)
	# large = (20, 100)
	# for comm_sizes in [large]:
	# 	for om in range(1, 6):
	# 		name = 'st_om_' + str(om)
	# 		LFR_graph_args[name] = {  # original
	# 			'N': 5000, 'k': 10, 'maxk': 50, 'minc': comm_sizes[0], 'maxc': comm_sizes[1],
	# 			't1': 2, 't2': 1, 'mu': 0.3, 'on': 2500, 'om': om}
	# 		LFR_graph_args[name] = {
	# 			'N': 1000, 'k': 10, 'maxk': 50, 'minc': comm_sizes[0], 'maxc': comm_sizes[1],
	# 			't1': 2, 't2': 1, 'mu': 0.3, 'on': 500, 'om': om}
	#
	# # Scale overlap
	# for overlap in range(0, 101, 10):
	# 	N = 2000
	# 	on = N / 100 * overlap
	# 	LFR_graph_args['on_' + str(overlap/100)[:3].rjust(3, '0')] = {
	# 		'N': 2000, 'k': 20, 'maxk': 50, 'minc': 60, 'maxc': 100,
	# 		't1': 2, 't2': 1, 'mu': 0.25, 'on': on, 'om': 2}
	#
	# # Scale mixing factor
	# for mu_factor in range(40, 71, 10):
	# 	om = 2
	# 	mu = 0.01 * mu_factor
	# 	k = om * 15  # Number of neighbors per community independent of mixing factor
	# 	k /= (1 - mu)
	# 	maxk = 1.5 * k  # Scale max degree with average degree
	# 	LFR_graph_args['mu_' + str(mu_factor).rjust(2,'0')] = {
	# 		'N': 2000, 'k': k, 'maxk': maxk, 'minc': 60, 'maxc': 100,
	# 		't1': 2, 't2': 2, 'mu': 0.01 * mu_factor, 'on': 2000, 'om': om}

	create_LFR_graphs(graphs, LFR_graph_args, iterations)
	for graph in graphs:
		graph.graph.indexEdges()
	return graphs


# ****************************************************************************************
# *                                     Algorithms                                       *
# ****************************************************************************************
def get_algos(storeEgoNets):
	algos = []
	algos.append(GroundTruth())
	# for iterations in [20, 40, 60, 80, 100]:
	# 	for threshold in [0.1]:
	# 		algos.append(SlpaAlgorithm(name='SLPA_{}_{}'.format(threshold, iterations),
	# 		                           threshold=threshold, numIterations=iterations))
	# algos.append(OlpAlgorithm())
	# algos.append(GceAlgorithm(alpha=1.0, add_name='_1.0'))
	# algos.append(GceAlgorithm(alpha=1.5, add_name='_1.5'))
	# algos.append(GceAlgorithm(alpha=1.0, add_name='_1.0_clean', clean_up='clean-merge'))
	# algos.append(MosesAlgorithm())
	# algos.append(OslomAlgorithm())

	partition_algos = OrderedDict()
	# partition_algos['PLP'] = [lambda g: PLP(g, 1, 20).run().getPartition()]
	# partition_algos['PLM'] = [lambda g: PLM(g, True, 1.0, 'none').run().getPartition()]
	# partition_algos['OSLOM'] = [lambda g: convertCoverToPartition(g, clusterOSLOM(g))]

	# def PLM_OSLOM_clean(g):
	# 	part = PLM(g, True, 1.0, 'none').run().getPartition()
	# 	cover_clean = OslomCleanUp(g, Cover(part)).run().getCover()
	# 	return convertCoverToPartition(g, cover_clean)
	# partition_algos['PLM_clean'] = [PLM_OSLOM_clean]
	# partition_algos['SLPA'] = [lambda g: SLPA(g, 0.1, 100).run().getPartition()]
	# partition_algos['PLM_merge'] = [
	# 	lambda g: partition_merge(g, PLM(copy(g), True, 1.0, 'none').run().getPartition())]
	# partition_algos['PLM-PLP_Info'] = [
	# 	lambda g: PLP(g, baseClustering=PLM(copy(g), True, 1.0, 'none').run().getPartition(),
	# 	              updateThreshold=g.numberOfNodes()/20, maxIterations=20).run().getPartition(),
	# 	lambda g: clusterInfomap(g)]
	# partition_algos['Potts_0.01'] = [lambda g: LPPotts(g, 0.01, 1, 20).run().getPartition()]
	# partition_algos['Potts_0.05'] = [lambda g: LPPotts(g, 0.05, 1, 20).run().getPartition()]
	# partition_algos['Potts_0.1'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition()]
	# partition_algos['LPPotts_par'] = [
	# 	lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition(),
	# 	lambda g: LPPotts(g, 0, 1, 20, True).run().getPartition()]
	# partition_algos['Infomap'] = [lambda g: clusterInfomap(g)]
	# partition_algos['Surprise'] = [lambda g: partitionLeiden(g, 'surprise')]
	partition_algos['Leiden_Mod'] = [lambda g: partitionLeiden(g, 'modularity')]
	# partition_algos['Leiden_Sig'] = [lambda g: partitionLeiden(g, 'significance')]
	# partition_algos['Leiden_SigRes'] = [lambda g: leidenSignificance(g)]

	for p_algos in partition_algos.values():
		p_algos.insert(1, lambda g: clusterInfomap(g))
	ego_parameters = get_ego_parameters(storeEgoNets)

	# TODO: Cleanup der overlapping groups verwirft statt merged
	algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='no-clean')
	algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='merge-overl')
	algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='remv-overl')
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='clean-merge')
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='clean-merge,overl')
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='clean_full')
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='clean-remove')
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='clean-merge-3')
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='clean-merge-5')
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='clean-merge-shrink')
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='clean-check')
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='clean-keep')
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='clean-keep-merge_E')
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='trim_gt')
	# algos += create_egosplit_algorithms(partition_algos, ego_parameters, clean_up='trim_gt-merge_gt')

	return algos


def get_ego_parameters(store_ego_nets):
	ego_parameters = OrderedDict()
	standard = {
		'weightFactor': 0,
		'weightOffset': 1,
		'storeEgoNet': 'Yes' if store_ego_nets else 'No',
		'addEgoNode': 'No',
		'extendStrategy': 'none',
		'extendPartitionIterations': 1,
		'addNodesFactor': 0,
		'addNodesExponent': 0,
		'partitionFromGroundTruth': 'No',
		'connectPersonas': 'Yes',
		'normalizePersonaCut': 'No',
		'connectPersonasStrat': 'spanning',
		'maxPersonaEdges': 1,
		'personaEdgeWeightFactor': 1,
		'normalizePersonaWeights': 'unweighted',
		'iterationWeight': 'No',
	}
	extend_standard = {
		**standard,
		'extendPartitionIterations': 1,
		'addNodesFactor': 2,
		'addNodesExponent': 0.8,
		'edgesBetweenNeigNeig': 'Yes',
		'minNodeDegree': 2,
		'extendOverDirected': 'No',
		'keepOnlyTriangles': 'No',
		'scoreStrategy': 'score',
		'onlyDirectedCandidates': 'No',
		'extendDirectedBack': 'Yes',
	}
	edge_scores_standard = {
		**extend_standard,
		'extendStrategy': 'edgeScore',
	}
	significance_scores_standard = {
		**edge_scores_standard,
		'extendPartitionIterations': 2,
		'extendStrategySecond': 'significance',
		'maxSignificance': 0.1,
		'sortGroups': 'significance',
		'orderedStatPos': 0.0,
		'subtractNodeDegree': 'Yes',
		'maxGroupsConsider': 10,
		'signMerge': 'Yes',
		'useSigMemo': 'No',
		'minEdgesToGroupSig': '3',  # TODO: should be 1
		'sigSecondRoundStrat': 'updateCandidates',
		'secondarySigExtRounds': '2',
	}

	# ego_parameters['b!'] = standard
	# ego_parameters['gt'] = {
	# 	**standard,
	# 	'partitionFromGroundTruth': 'Yes',
	# }
	# ego_parameters['e-2-0.8'] = {
	# 	**edge_scores_standard,
	# 	'addNodesFactor': 2,
	# }
	ego_parameters['e!'] = {
		**edge_scores_standard,
	}
	# ego_parameters['b-s*1'] = {
	# 	**significance_scores_standard,
	# 	'extendStrategy': 'none',
	# }
	# ego_parameters['b-s*3'] = {
	# 	**significance_scores_standard,
	# 	'extendStrategy': 'none',
	# 	'extendPartitionIterations': 4,
	# }
	# ego_parameters['e-s*1'] = {
	# 	**significance_scores_standard,
	# }
	# ego_parameters['e-s*3'] = {
	# 	**significance_scores_standard,
	# 	'extendPartitionIterations': 4,
	# }
	# ego_parameters['s-order'] = {
	# 	**significance_scores_standard,
	# 	'useSigMemo': 'No',
	# 	'sigSecondRoundStrat': 'orderStat',
	# }
	# ego_parameters['s-mem-update'] = {
	# 	**significance_scores_standard,
	# 	'useSigMemo': 'Yes',
	# 	'sigSecondRoundStrat': 'updateCandidates',
	# }
	# ego_parameters['s-mem'] = {
	# 	**significance_scores_standard,
	# 	'useSigMemo': 'Yes',
	# }
	# ego_parameters['s-2'] = {
	# 	**significance_scores_standard,
	# 	'addNodesFactor': 1,
	# 	'orderedStatPos': 0.1,
	# 	'useSigMemo': 'No',
	# 	'minEdgesToGroupSig': '1',
	# }

	return ego_parameters


def create_LFR_graphs(graphs, graphArgs, iterations):
	for argsName in graphArgs.keys():
		args = graphArgs[argsName]
		graphName = 'LFR_' + argsName
		for i in range(iterations):
			graph_wrapper = LFRGraph(graphName, args)
			graphs.append(graph_wrapper)
	return graphs


def create_egosplit_algorithms(partition_algos, ego_parameters, clean_up=''):
	algos = []
	for part_name in partition_algos:
		for para_name, parameters in ego_parameters.items():
			name = '{}{}{}{}'.format('Ego_', part_name,
			                         '_' + para_name if para_name else '',
			                         '_' + clean_up if clean_up else '')
			algos.append(EgoSplitAlgorithm(
				name,
				{key.encode('utf-8'): str(value).encode('utf-8')
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
