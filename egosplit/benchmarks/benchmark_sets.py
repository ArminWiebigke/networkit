from egosplit.benchmarks.data_structures.benchmark_set import BenchmarkSet
from egosplit.benchmarks.execution.algo_creation import	OtherAlgorithms
from egosplit.benchmarks.execution.egosplit_config import EgoSplitClusteringAlgorithmsConfig, EgoSplitParameterConfig
from egosplit.benchmarks.execution.cleanup import CleanUpConfig
from egosplit.benchmarks.execution.graph_creation import GraphSetsConfig
from egosplit.benchmarks.plot_scripts.bench_config import PlotGraphSetConfig, PlotAlgoSetConfig
from egosplit.benchmarks.plot_scripts.create_plots import PlotSetConfig


def get_benchmark_configs():
	benchmark_sets = [
		# Scratchpad(),
		# EdgesScore(),
		# EdgesScoreSignificance(),
		# EdgesFactor(),
		# SigMerge(),
		# SigExtIter(),
		# SigCheckUpdated(),
		# SigMaxCandidates(),
		# SigClusterIter(),
		# SigMemoize(),
		# ExtensionCompare(),
		# LocalClustering(),
		# LocalClusteringPlots(),
		# ConnectPersonas(),
		# GlobalClustering(),
		# CleanUp(),
		# NewCleanUp(),
		# CleanUpSnap(),
		# CompareOther(),
		# EgoSnap(),

		LocalClusteringBest(),
		GlobalClusteringBest(),
		OptimizationSteps(),
		AlgorithmsSnap(),
		AlgorithmsLFR(),
	]
	# TODO: Wonanders hin?
	# To enforce a specific order of the algorithms, numbers are added to algorithm names.
	# In the plots, these numbers should be removed.
	for benchmark_set in benchmark_sets:
		benchmark_set.remove_algo_parts.extend(['({:03.0f})'.format(i) for i in range(30)])
	return benchmark_sets


class Scratchpad(BenchmarkSet):
	config = {
		'name':
			'test',
		'result_dir':
			'test',
		'plot_dir':
			'test/',
		'time_limit':
			1000,
		'iterations': 1,
		EgoSplitClusteringAlgorithmsConfig:
			'test',
		EgoSplitParameterConfig:
			['test'],
		CleanUpConfig: 'No Cleanup',
		# 'stream_to_gephi':
		# 	True,
		# 'store_ego_nets':
		# 	True,
		GraphSetsConfig: [
			# 'om',
			# 'mu',
			# 'facebook',
			# 'overlap',
			# large_graphs,
			'test',
			# 'snap',
		],
		PlotGraphSetConfig: [
			# 'om',
			# 'facebook_bar',
			# 'facebook',
		],
		PlotSetConfig: [
			'metrics',
			'timings',
			# 'ego_net_cluster',
			# 'comm_sizes',
			# 'comm_f1',
		],
		PlotAlgoSetConfig: [
			'all'
		],
		'remove_algo_parts':
			['Ego', ' + Infomap', ' | ', 'No Clean Up', 'No Extension',
			 'EdgesScore', 'Significance'],
	}


class LocalClusteringBest(BenchmarkSet):
	config = {
		'name':
			'local_cluster_best',
		'result_dir':
			'local_cluster_best',
		'plot_dir':
			'local_cluster_best/',
		EgoSplitClusteringAlgorithmsConfig:
			'local',
		EgoSplitParameterConfig:
			['edges'],
		CleanUpConfig:
			'best',
		GraphSetsConfig: [
			'om',
			'overlap',
			'mu',
			# 'facebook',
			# 'snap',
		],
		PlotGraphSetConfig: [
			'om',
			'mu',
			# 'facebook',
			# 'facebook_bar',
			# 'snap',
		],
		PlotSetConfig: [
			'metrics',
		],
		PlotAlgoSetConfig: [
			'all'
		],
		'remove_algo_parts':
			['Ego', ' | ', ' + LM-Map', 'Clean-new', 'EdgesScore'],
		'replace_legend':
			{'Leiden': 'LeidenMod', 'Potts': 'LPPotts'},
	}


class GlobalClusteringBest(BenchmarkSet):
	config = {
		'name':
			'global_cluster_best',
		'result_dir':
			'global_cluster_best',
		'plot_dir':
			'global_cluster_best/',
		EgoSplitClusteringAlgorithmsConfig:
			'global',
		EgoSplitParameterConfig:
			['edges'],
		CleanUpConfig:
			'best',
		'calc_f1_per_comm':
			True,
		GraphSetsConfig: [
			'om',
			'overlap',
			'mu',
			# 'facebook',
			# 'snap',
			# 'test',
		],
		PlotGraphSetConfig: [
			'om',
			'mu',
			# 'facebook',
			# 'facebook_bar',
			# 'snap',
		],
		PlotSetConfig: [
			'metrics',
			'timings',
			'ego_net_cluster',
			'comm_sizes',
			'comm_f1',
		],
		PlotAlgoSetConfig: [
			'all'
		],
		'remove_algo_parts':
			['Ego', ' | ', 'Clean-new', 'EdgesScore', 'Leiden + '],
		'replace_legend':
			{'Leiden': 'LeidenMod', 'Potts': 'LPPotts'},
	}


class OptimizationSteps(BenchmarkSet):
	config = {
		'name':
			'optimization_steps',
		'result_dir':
			'optimization_steps',
		'plot_dir':
			'optimization_steps/',
		EgoSplitClusteringAlgorithmsConfig:
			'best',
		EgoSplitParameterConfig:
			['steps'],
		GraphSetsConfig: [
			'om',
			'overlap',
			'mu',
		],
		PlotGraphSetConfig: [
			'om',
			'mu',
		],
		PlotSetConfig: [
			'metrics',
		],
		PlotAlgoSetConfig: [
			'all'
		],
		'remove_algo_parts':
			['Ego', ' | ', 'Leiden + LM-Map'],
		'replace_legend':
			{'Base': 'Ego', 'EdgesScore': 'Ego+E', 'E + Spanning': 'Ego+E+S', 'E + S + Cleanup': 'Ego+E+S+C'},
	}


class EgoSnap(BenchmarkSet):
	config = {
		'name':
			'ego_snap',
		'result_dir':
			'ego_snap',
		'plot_dir':
			'ego_snap/',
		'iterations': 1,
		'time_limit': 7200,
		EgoSplitClusteringAlgorithmsConfig:
			'two_best',
		EgoSplitParameterConfig:
			['best'],
		GraphSetsConfig: [
			'snap',
		],
		OtherAlgorithms: [
			'Ego-original',
		],
		PlotGraphSetConfig: [
			'snap'
		],
		PlotSetConfig: [
			'metrics',
			'timings',
		],
		PlotAlgoSetConfig:
			['all'],
		'remove_algo_parts':
			['Ego', ' |', 'EdgesScore'],
		'replace_legend':
			{'No Clean Up': 'Uncleaned', 'Clean-new': 'Cleaned'},
	}


class AlgorithmsLFR(BenchmarkSet):
	config = {
		'name':
			'other_algos_lfr',
		'result_dir':
			'other_algos_lfr',
		'plot_dir':
			'other_algos_lfr/',
		'time_limit':
			1800,
		OtherAlgorithms: [
			'Ego-original',
			'Ego-optimized',
			'Oslom',
			'Moses',
			'GCE',
			'Demon'
		],
		GraphSetsConfig: [
			'om',
			'overlap',
			'mu',
		],
		PlotGraphSetConfig: [
			'om',
			'mu',
		],
		PlotSetConfig: [
			'metrics',
			'timings',
		],
		PlotAlgoSetConfig: [
			'all'
		],
		'replace_legend':
			{'Ego-original': 'Ego', 'Ego-optimized': 'Ego+'},
	}


class EgoSplittingLFR(BenchmarkSet):
	config = {
		'name':
			'egosplitting_lfr',
		'result_dir':
			'egosplitting_lfr',
		'plot_dir':
			'egosplitting_lfr/',
		'iterations':
			10,
		'time_limit':
			3600,
		EgoSplitClusteringAlgorithmsConfig:
			'plm_lm-map',
		EgoSplitParameterConfig:
			['CleanupMinOverlap'],
		OtherAlgorithms: [
			'Ego-original',
		],
		GraphSetsConfig: [
			'om',
                        'overlap',
                        'mu',
                        #'test',
		],
		PlotGraphSetConfig: [
			'om',
			'mu',
		],
		PlotSetConfig: [
			'metrics',
			'timings',
		],
		PlotAlgoSetConfig: [
			'all'
		],
		'remove_algo_parts':
			['Ego', ' | ', 'PLM + LM-Map'],
	}



class EgoSplittingSnap(BenchmarkSet):
	config = {
		'name':
			'egosplitting_snap',
		'result_dir':
			'egosplitting_snap',
		'plot_dir':
			'egosplitting_snap/',
		'iterations':
			10,
		'time_limit':
			7200,
		EgoSplitClusteringAlgorithmsConfig:
			'plm_lm-map',
		EgoSplitParameterConfig:
			['CleanupMinOverlap'],
		OtherAlgorithms: [
			'Ego-original',
		],
		GraphSetsConfig: [
			'snap',
                        #'test',
		],
		PlotGraphSetConfig: [
			'snap',
		],
		PlotSetConfig: [
			'metrics',
			'timings',
		],
		PlotAlgoSetConfig: [
			'all'
		],
		'remove_algo_parts':
			['Ego', ' | ', 'PLM + LM-Map'],
		'replace_legend':
			{'Ego-original': 'Ego', 'Ego-optimized': 'Ego+'},
	}



class AlgorithmsSnap(BenchmarkSet):
	config = {
		'name':
			'other_algos_snap',
		'result_dir':
			'other_algos_snap',
		'plot_dir':
			'other_algos_snap/',
		'iterations':
			1,
		'time_limit':
			7200,
		OtherAlgorithms: [
			'Ego-original',
			'Ego-optimized',
			'Oslom',
			'Moses',
			'GCE',
			'Demon'
		],
		GraphSetsConfig: [
			'snap',
		],
		PlotGraphSetConfig: [
			'snap',
		],
		PlotSetConfig: [
			'metrics',
			'timings',
		],
		PlotAlgoSetConfig: [
			'all'
		],
		'replace_legend':
			{'Ego-original': 'Ego', 'Ego-optimized': 'Ego+'},
	}


class EdgesScore(BenchmarkSet):
	config = {
		'name':
			'edges_score',
		'result_dir':
			'edges_score',
		'plot_dir':
			'extend/edges/score/',
		EgoSplitClusteringAlgorithmsConfig:
			'fast',
		EgoSplitParameterConfig:
			['no-extend', 'edges-score'],
		'store_ego_nets':
			True,
		'iterations':
			1,
		GraphSetsConfig: [
			'om',
		],
		PlotGraphSetConfig: [
			'om',
		],
		PlotSetConfig: [
			'ego_net_extend',
			'timings',
			'metrics'
		],
		PlotAlgoSetConfig: [
			'all',
		],
		'remove_algo_parts':
			['Ego', ' | ', 'PLP + PLM', 'No Clean Up'],
		'replace_legend': {
			# 'Extend: Edges': '$q_1$',
			# 'Extend: Edges div Degree': '$q_2$',
			# 'Extend: Edges pow 2 div Degree': '$q_3$',
			# 'Extend: Random': '$q_4$'
		},
	}


class EdgesScoreSignificance(BenchmarkSet):
	config = {
		'name':
			'edges_score_sig',
		'result_dir':
			'edges_score_sig',
		'plot_dir':
			'extend_edges_significance/',
		EgoSplitClusteringAlgorithmsConfig:
			'best',
		EgoSplitParameterConfig:
			['no-extend', 'edges-significance'],
		'store_ego_nets':
			True,
		GraphSetsConfig: [
			'om',
			# 'mu',
			'facebook'
		],
		PlotGraphSetConfig: [
			'om',
			# 'mu',
			'facebook',
		],
		PlotSetConfig: [
			'ego_net_extend',
			'timings',
			'metrics'
		],
		PlotAlgoSetConfig: [
			'Leiden-Info',
		],
		'remove_algo_parts':
			['Ego', ' | ', 'Leiden + Infomap', 'No Clean Up'],
	}


class EdgesFactor(BenchmarkSet):
	config = {
		'name':
			'edges_factor',
		'result_dir':
			'edges_factor',
		'plot_dir':
			'extend/edges/add_factor/',
		EgoSplitClusteringAlgorithmsConfig:
			'fast',
		EgoSplitParameterConfig:
			['no-extend', 'edges-factor'],
		'store_ego_nets':
			True,
		'score_per_egonet':
			True,
		GraphSetsConfig:
			['om'],
		PlotGraphSetConfig:
			['om'],
		PlotSetConfig:
			['ego_net_extend', 'ego_net_x_extend', 'timings'],
		PlotAlgoSetConfig:
			['all'],
		'remove_algo_parts':
			['Ego', ' | ', 'PLP + PLM', 'No Clean Up'],
	}


class SigMerge(BenchmarkSet):
	config = {
		'name':
			'sig_merge',
		'result_dir':
			'sig_merge',
		'plot_dir':
			'extend/sig/merge_groups/',
		EgoSplitClusteringAlgorithmsConfig:
			'leiden local',
		EgoSplitParameterConfig:
			['no-extend', 'sig-merge'],
		'store_ego_nets':
			True,
		GraphSetsConfig:
			['om'],
		PlotGraphSetConfig:
			['om'],
		PlotSetConfig:
			['ego_net_extend', 'timings'],
		PlotAlgoSetConfig:
			['all'],
		'remove_algo_parts':
			['Ego', ' | ', 'Leiden + Infomap', 'No Clean Up'],
		'replace_legend':
			{'Single Clusters': 'Single', 'Single + Merged Clusters': 'Merged'},
	}


class SigExtIter(BenchmarkSet):
	config = {
		'name':
			'sig_ext_iter',
		'result_dir':
			'sig_ext_iter',
		'plot_dir':
			'extend/sig/extend_iterative/',
		EgoSplitClusteringAlgorithmsConfig:
			'leiden local',
		EgoSplitParameterConfig:
			['no-extend', 'sig-ext-iter'],
		'store_ego_nets':
			True,
		GraphSetsConfig:
			['om'],
		PlotGraphSetConfig:
			['om'],
		PlotSetConfig:
			['ego_net_extend', 'timings'],
		PlotAlgoSetConfig:
			['all'],
		'remove_algo_parts':
			['Ego', ' | ', 'Leiden + Infomap', 'No Clean Up', '!'],
		'replace_legend': {
			'0 Iterations': '$I_{max} = 0$',
			'1 Iteration': '$I_{max} = 1$',
			'2 Iterations': '$I_{max} = 2$',
			'3 Iterations': '$I_{max} = 3$',
			'5 Iterations': '$I_{max} = 5$',
			'10 Iterations': '$I_{max} = 10$',
			'100 Iterations': '$I_{max} = 100$',
		},
	}


class SigCheckUpdated(BenchmarkSet):
	config = {
		'name':
			'sig_check_updated',
		'result_dir':
			'sig_check_updated',
		'plot_dir':
			'extend/sig/check_updated/',
		EgoSplitClusteringAlgorithmsConfig:
			'leiden local',
		EgoSplitParameterConfig:
			['no-extend', 'sig-check-updated'],
		'store_ego_nets':
			True,
		GraphSetsConfig:
			['om'],
		PlotGraphSetConfig:
			['om'],
		PlotSetConfig:
			['ego_net_extend', 'timings'],
		PlotAlgoSetConfig:
			['all'],
		'remove_algo_parts':
			['Ego', ' | ', 'Leiden + Infomap', 'No Clean Up'],
	}


class SigMaxCandidates(BenchmarkSet):
	config = {
		'name':
			'sig_max_candidates',
		'result_dir':
			'sig_max_candidates',
		'plot_dir':
			'extend/sig/max_candidates/',
		EgoSplitClusteringAlgorithmsConfig:
			'leiden local',
		EgoSplitParameterConfig:
			['no-extend', 'sig-max-candidates'],
		'store_ego_nets':
			True,
		GraphSetsConfig:
			['om'],
		PlotGraphSetConfig:
			['om'],
		PlotSetConfig:
			['ego_net_extend', 'timings'],
		PlotAlgoSetConfig:
			['all'],
		'remove_algo_parts':
			['Ego', ' | ', 'Leiden + Infomap', 'No Clean Up'],
	}


class SigClusterIter(BenchmarkSet):
	config = {
		'name': 'sig_cluster_iter',
		'result_dir': 'sig-cluser-iter',
		'plot_dir': 'extend/sig/cluster_iterative/',
		EgoSplitClusteringAlgorithmsConfig: 'leiden local',
		EgoSplitParameterConfig: ['no-extend', 'sig-cluster-iter'],
		'store_ego_nets': True,
		GraphSetsConfig:
			['om'],
		PlotGraphSetConfig:
			['om'],
		PlotSetConfig:
			['ego_net_extend', 'timings'],
		PlotAlgoSetConfig:
			['all'],
		'remove_algo_parts': ['Ego', ' | ', 'Leiden + Infomap', 'No Clean Up'],
	}


class SigMemoize(BenchmarkSet):
	config = {
		'name': 'sig_mem',
		'result_dir': 'sig-mem',
		'plot_dir': 'extend/sig/memoize/',
		EgoSplitClusteringAlgorithmsConfig: 'leiden local',
		EgoSplitParameterConfig: ['no-extend', 'sig-mem'],
		'store_ego_nets': True,
		GraphSetsConfig:
			['om'],
		PlotGraphSetConfig:
			['om'],
		PlotSetConfig:
			['ego_net_extend', 'timings'],
		PlotAlgoSetConfig:
			['all'],
		'remove_algo_parts': ['Ego', ' | ', 'Leiden + Infomap', 'No Clean Up'],
	}


class ExtensionCompare(BenchmarkSet):
	config = {
		'name':
			'ext_compare',
		'result_dir':
			'ext_compare',
		'plot_dir':
			'extend/compare/',
		EgoSplitClusteringAlgorithmsConfig:
			'leiden local',
		EgoSplitParameterConfig:
			['no-extend', 'extend'],
		'store_ego_nets':
			True,
		GraphSetsConfig: [
			'om',
			'overlap',
			# 'mu',
			'facebook',
			'snap',
		],
		PlotGraphSetConfig: [
			'om',
			# 'mu',
			'facebook',
			# 'facebook_bar',
			'snap',
		],
		PlotSetConfig: [
			'ego_net_extend',
			'timings',
		],
		PlotAlgoSetConfig:
			['all'],
		'remove_algo_parts':
			['Ego', ' | ', 'Leiden + Infomap', 'No Clean Up'],
	}


class LocalClustering(BenchmarkSet):
	config = {
		'name':
			'local_cluster',
		'result_dir':
			'local_cluster',
		'plot_dir':
			'local_cluster/',
		EgoSplitClusteringAlgorithmsConfig:
			'local',
		EgoSplitParameterConfig:
			['no-extend', 'extend'],
		'store_ego_nets':
			True,
		GraphSetsConfig: [
			'om',
			'overlap',
			# 'mu',
			'facebook',
			# 'snap',
		],
		PlotGraphSetConfig: [
			'om',
			# 'mu',
			'facebook',
			# 'facebook_bar',
			# 'snap',
		],
		PlotSetConfig: [
			'metrics',
			'timings',
			'ego_net_cluster',
			'comm_sizes',
			'comm_f1',
		],
		PlotAlgoSetConfig: [
			'base',
			'edges',
			'sig',
		],
		'remove_algo_parts':
			['Ego', ' + MapEquation', ' | ', 'No Clean Up', 'No Extension',
			 'EdgesScore', 'Significance'],
		'replace_legend':
			{'Leiden': 'LeidenMod', 'Potts': 'LPPotts'},
	}


class LocalClusteringPlots(BenchmarkSet):
	config = {
		'name':
			'local_cluster',
		'result_dir':
			'local_cluster',
		'plot_dir':
			'local_cluster_per_algo/',
		EgoSplitClusteringAlgorithmsConfig:
			'local',
		EgoSplitParameterConfig:
			['no-extend', 'extend'],
		'store_ego_nets':
			True,
		GraphSetsConfig: [
		],
		PlotGraphSetConfig: [
			# 'om',
			# 'mu',
			# 'facebook',
			# 'facebook_bar',
			'snap',
		],
		PlotSetConfig: [
			'metrics',
			'timings',
			'ego_net_extend',
			'ego_net_cluster',
			'comm_sizes',
			'comm_f1',
			'num_comms',
		],
		PlotAlgoSetConfig: [
			'PLP',
			'PLM',
			'Potts',
		],
		'remove_algo_parts':
			['Ego', ' + Infomap', ' | ', 'No Clean Up'],
		'replace_legend':
			{'Leiden': 'LeidenMod', 'Potts': 'LPPotts'},
	}


class ConnectPersonas(BenchmarkSet):
	config = {
		'name':
			'connect_persona',
		'result_dir':
			'connect_persona',
		'plot_dir':
			'connect_persona/',
		EgoSplitClusteringAlgorithmsConfig:
			'Leiden/Infomap + Infomap',
		EgoSplitParameterConfig:
			['connect-persona'],
		'store_ego_nets':
			True,
		GraphSetsConfig:
			['om', 'overlap', 'mu', 'facebook'],
		PlotGraphSetConfig: [
			'om',
			'mu',
			'facebook',
			'facebook_bar'
		],
		PlotSetConfig: [
			'metrics',
			'timings',
			'comm_f1',
			'comm_sizes',
		],
		PlotAlgoSetConfig:
			['connect-persona'],
		'remove_algo_parts':
			['Ego', ' | ', 'Leiden + Infomap', 'Infomap + Infomap', 'No Clean Up', 'No Extension',
			 'EdgesScore'],
		'replace_legend':
			{'No Connection': 'NoConnection', 'Max Spanning Unweighted': 'SpanUnweight',
			 'All Unweighted': 'AllUnweight', 'All Density Max Weight 1': 'AllWeight'},
	}


class GlobalClustering(BenchmarkSet):
	config = {
		'name':
			'global_cluster',
		'result_dir':
			'global_cluster',
		'plot_dir':
			'global_cluster/',
		EgoSplitClusteringAlgorithmsConfig:
			'global',
		EgoSplitParameterConfig:
			['edges'],
		GraphSetsConfig: [
			'om',
			'overlap',
			# 'mu',
			'facebook',
			'snap',
		],
		PlotGraphSetConfig: [
			'om',
			# 'mu',
			'facebook',
			# 'facebook_bar'
			'snap',
		],
		PlotSetConfig: [
			'metrics',
			'timings',
			'comm_f1',
			'comm_sizes',
		],
		PlotAlgoSetConfig:
			['Info-local', 'Leiden-local'],
		'remove_algo_parts':
			['Ego', ' | ', 'No Clean Up', 'Infomap + ', 'Leiden + ', 'No Extension', 'EdgesScore'],
		'replace_legend':
			{'Leiden': 'LeidenMod', 'Potts': 'LPPotts'},
	}


class CleanUp(BenchmarkSet):
	config = {
		'name':
			'clean_up',
		'result_dir':
			'clean_up',
		'plot_dir':
			'clean_up/',
		EgoSplitClusteringAlgorithmsConfig:
			'two_best',
		EgoSplitParameterConfig:
			['edges'],
		CleanUpConfig:
			'all',
		GraphSetsConfig: [
			'om',
			'overlap',
			'mu',
			'facebook',
		],
		PlotGraphSetConfig: [
			'om',
			'mu',
			'facebook',
			'facebook_bar'
		],
		PlotSetConfig: [
			'metrics',
			'timings',
			'comm_f1',
			'comm_sizes',
			'num_comms',
		],
		PlotAlgoSetConfig:
			['clean-up'],
		'remove_algo_parts':
			['Ego', ' | ', 'EdgesScore', 'Infomap + Surprise',
			 'Leiden + Infomap',
			 ],
		'replace_legend':
			{'No Clean Up': 'Original', 'Clean-merge': 'Cleaned',
			 'Clean-remove': 'CleanRemove', 'OSLOM-full': 'OSLOM'},
	}


class NewCleanUp(BenchmarkSet):
	config = {
		'name':
			'new_clean_up',
		'result_dir':
			'new_clean_up',
		'plot_dir':
			'new_clean_up/',
		EgoSplitClusteringAlgorithmsConfig:
			'best',
		EgoSplitParameterConfig:
			['edges'],
		CleanUpConfig:
			'new_clean',
		GraphSetsConfig: [
			'om',
			'overlap',
			# 'mu',
			# 'facebook',
		],
		PlotGraphSetConfig: [
			'om',
			# 'mu',
			# 'facebook',
			# 'facebook_bar'
		],
		PlotSetConfig: [
			'metrics',
			'timings',
			'comm_f1',
			'comm_sizes',
			'num_comms',
		],
		PlotAlgoSetConfig:
			['all'],
		'remove_algo_parts':
			['Ego', ' | ', 'EdgesScore', 'Infomap + Surprise',
			 'Leiden + Infomap',
			 ],
		'replace_legend':
			{'No Clean Up': 'Uncleaned', 'Clean-merge': 'Cleaned old',
			 'Clean-new': 'Cleaned new'},
	}


class CleanUpSnap(BenchmarkSet):
	config = {
		'name':
			'clean_up_snap',
		'result_dir':
			'clean_up_snap',
		'plot_dir':
			'clean_up_snap/',
		EgoSplitClusteringAlgorithmsConfig:
			'two_best',
		EgoSplitParameterConfig:
			['edges'],
		'other_algos': [
			'Ego-original',
		],
		CleanUpConfig:
			'with-without',
		GraphSetsConfig: [
			'snap',
		],
		PlotGraphSetConfig: [
			'snap'
		],
		PlotSetConfig: [
			'metrics',
			'timings',
			'comm_f1',
			'comm_sizes',
			'num_comms',
		],
		PlotAlgoSetConfig:
			['all'],
		'remove_algo_parts':
			['Ego', ' |', 'EdgesScore'],
		'replace_legend':
			{'No Clean Up': 'Uncleaned', 'Clean-new': 'Cleaned'},
	}


class CompareOther(BenchmarkSet):
	config = {
		'name':
			'compare_other',
		'result_dir':
			'compare_other',
		'plot_dir':
			'compare_other/',
		EgoSplitClusteringAlgorithmsConfig:
			'best',
		EgoSplitParameterConfig:
			['edges'],
		CleanUpConfig:
			'two_best',
		'other_algos': [
			'GCE',
			'Moses',
			'Oslom',
			'Ego-original',
			'Ego-base',
		],
		GraphSetsConfig: [
			# 'om',
			'overlap',
			# 'mu',
			# 'facebook',
			# 'snap',
		],
		PlotGraphSetConfig: [
			'om',
			# graph_sets['om_max'],
			# 'mu',
			# graph_sets['mu_max'],
			# 'facebook',
			# 'facebook_bar'
		],
		PlotSetConfig: [
			'metrics',
			'timings',
		],
		PlotAlgoSetConfig:
			['all'],
		'remove_algo_parts':
			[' | ', 'Leiden + Infomap', 'EdgesScore', 'Clean-merge'],
		'replace_legend':
			{'Ego': 'ES+', 'Ego-original': 'ESF', 'Ego-base': 'ESB'},
	}
