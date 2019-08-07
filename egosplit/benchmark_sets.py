
def get_benchmark_configs():
	for config in benchmark_configs.values():
		set_default_args(config)

	configs = []
	for b in benchmark_sets:
		configs.append(benchmark_configs[b])
	return configs


benchmark_sets = [
	# 'test',
	# 'edges-score',
	# 'edges-factor',
	# 'sig-merge',
	# 'sig-max-candidates',
	# 'sig-ext-iter',
	# 'sig-check-updated',
	# 'sig-cluster-iter',
	# 'ext-compare',
	# 'local-cluster',
	# 'connect-persona',
	# 'global-cluster',
	'clean-up',
]


def set_default_args(config):
	default_args = {
		'clean_up_set': 'No Cleanup',
		'hue': 'Algorithm',
		'score_per_egonet': False,
		'store_ego_nets': False,
		'other_algos': None,
	}
	for key, value in default_args.items():
		if key not in config:
			config[key] = value


benchmark_configs = {
	'edges-score': {
		'result_dir': 'edges-score',
		'plot_dir': 'extend/edges/score_norm/',
		'ego_part_algos': 'fast',
		'ego_params': ['no-extend', 'edges-score'],
		'store_ego_nets': True,
		'graph_sets': [
			'om',
			# 'mu',
			# 'facebook'
		],
		'remove_algo_parts': ['Ego', 'PLP + PLM | ', ' | No Clean Up', '!',
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {'Extend: Edges': '$q_1$', 'Extend: Edges / Degree': '$q_2$',
		                   'Extend: Edges^2 / Degree': '$q_3$', 'Extend: Random': '$q_4$'},
		'plots': ['ego_net_extend'],
		'plot_algo_set': ['all'],

	},
	'edges-factor': {
		'result_dir': 'edges-factor',
		'plot_dir': 'extend/edges/add_factor/',
		'ego_part_algos': 'fast',
		'ego_params': ['no-extend', 'edges-factor'],
		'store_ego_nets': True,
		'score_per_egonet': True,
		'graph_sets': [
			'om',
			# 'mu',
			# 'facebook'
		],
		'remove_algo_parts': ['Ego', 'PLP + PLM | ', ' | No Clean Up', '!',
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {},
		'plots': ['ego_net_extend', 'ego_net_x_extend', 'timings'],
		'plot_algo_set': ['all'],
	},
	'sig-merge': {
		'result_dir': 'sig-merge',
		'plot_dir': 'extend/sig/merge_groups/',
		'ego_part_algos': 'leiden local',
		'ego_params': ['no-extend', 'sig-merge'],
		'store_ego_nets': True,
		'graph_sets': [
			'om',
			# 'mu',
			# 'facebook'
		],
		'remove_algo_parts': ['Ego', 'Leiden + Infomap | ', ' | No Clean Up', '!',
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {'Single Clusters': 'Single', 'Single + Merged Clusters': 'Merged'},
		'plots': ['ego_net_extend', 'timings'],
		'plot_algo_set': ['all'],
	},
	'sig-ext-iter': {
		'result_dir': 'sig-ext-iter',
		'plot_dir': 'extend/sig/extend_iterative/',
		'ego_part_algos': 'leiden local',
		'ego_params': ['no-extend', 'sig-ext-iter'],
		'store_ego_nets': True,
		'graph_sets': [
			'om',
			# 'mu',
			# 'facebook'
		],
		'remove_algo_parts': ['Ego', 'Leiden + Infomap | ', ' | No Clean Up', '!',
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {},
		'plots': ['ego_net_extend', 'timings'],
		'plot_algo_set': ['all'],
	},
	'sig-check-updated': {
		'result_dir': 'sig-check-updated',
		'plot_dir': 'extend/sig/check_updated/',
		'ego_part_algos': 'leiden local',
		'ego_params': ['no-extend', 'sig-check-updated'],
		'store_ego_nets': True,
		'graph_sets': [
			'om',
			# 'mu',
			# 'facebook'
		],
		'remove_algo_parts': ['Ego', 'Leiden + Infomap | ', ' | No Clean Up', '!',
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {},
		'plots': ['ego_net_extend', 'timings'],
		'plot_algo_set': ['all'],
	},
	'sig-max-candidates': {
		'result_dir': 'sig-max-candidates',
		'plot_dir': 'extend/sig/max_candidates/',
		'ego_part_algos': 'leiden local',
		'ego_params': ['no-extend', 'sig-max-candidates'],
		'store_ego_nets': True,
		'graph_sets': [
			'om',
			# 'mu',
			# 'facebook'
		],
		'remove_algo_parts': ['Ego', 'Leiden + Infomap | ', ' | No Clean Up', '!',
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {},
		'plots': ['ego_net_extend', 'timings'],
		'plot_algo_set': ['all'],
	},
	'sig-cluster-iter': {
		'result_dir': 'sig-cluser-iter',
		'plot_dir': 'extend/sig/cluster_iterative/',
		'ego_part_algos': 'leiden local',
		'ego_params': ['no-extend', 'sig-cluster-iter'],
		'store_ego_nets': True,
		'graph_sets': [
			'om',
			# 'mu',
			# 'facebook'
		],
		'remove_algo_parts': ['Ego', 'Leiden + Infomap | ', ' | No Clean Up', '!',
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {},
		'plots': ['ego_net_extend', 'timings'],
		'plot_algo_set': ['all'],
	},
	'ext-compare': {
		'result_dir': 'ext-compare',
		'plot_dir': 'extend/compare/',
		'ego_part_algos': 'leiden local',
		'ego_params': ['no-extend', 'extend'],
		'store_ego_nets': True,
		'score_per_egonet': True,
		'graph_sets': [
			'om',
			'mu',
			'facebook'
		],
		'remove_algo_parts': ['Ego', 'Leiden + Infomap | ', ' | No Clean Up', '!',
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {},
		# 'plots': ['ego_net_extend', 'ego_net_x_extend', 'timings'],
		'plots': ['ego_net_extend', 'metrics', 'timings'],
		'plot_algo_set': ['all'],
	},
	'local-cluster': {
		'result_dir': 'local-cluster',
		'plot_dir': 'local_cluster/',
		'ego_part_algos': 'local',
		'ego_params': ['no-extend', 'extend'],
		'store_ego_nets': True,
		'score_per_egonet': True,
		'graph_sets': [
			'om',
			'overlap',
			'mu',
			'facebook',
			'facebook_bar'
		],
		'remove_algo_parts': ['Ego', ' + Infomap', ' | No Clean Up', '!', ' | No Extension',
		                      ' | EdgesScore', ' | Significance'
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {},
		'plots': ['ego_net_cluster', 'metrics', 'timings'],
		'plot_algo_set': ['base', 'edges', 'sig'],
	},
	'connect-persona': {
		'result_dir': 'connect-persona',
		'plot_dir': 'connect_persona/',
		'ego_part_algos': "Leiden/Infomap + Infomap",
		'ego_params': ['connect-persona'],
		'store_ego_nets': True,
		'graph_sets': [
			'om',
			'overlap',
			'mu',
			'facebook',
		],
		'plot_graph_sets': [
			'om',
			'mu',
			# 'facebook',
			'facebook_bar',
		],
		'remove_algo_parts': ['Ego', 'Leiden + Infomap', 'Infomap + Infomap', ' | No Clean Up',
		                      ' | No Extension',
		                      ' | EdgesScore | ',
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {'No Connection': 'NoConnection',
		                   'Max Spanning Unweighted': 'MaxSpanUnweight',
		                   'All Unweighted': 'AllUnweight',
		                   'All Density Max Weight 1': 'AllWeight'},
		'plots': ['metrics', 'timings'],
		'plot_algo_set': ['Leiden-Info', 'Info-Info'],
	},
	'global-cluster': {
		'result_dir': 'global-cluster',
		'plot_dir': 'global_cluster/',
		'ego_part_algos': "global",
		'ego_params': ['edges'],
		'graph_sets': [
			'om',
			'overlap',
			'mu',
			# 'facebook',
		],
		'plot_graph_sets': [
			'om',
			'mu',
			'facebook',
			'facebook_bar',
		],
		'remove_algo_parts': ['Ego', ' | No Clean Up', 'Infomap + ', 'Leiden + ',
		                      ' | No Extension',
		                      ' | EdgesScore',
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {},
		'plots': ['metrics', 'comm_sizes', 'timings'],
		'plot_algo_set': ['Info-local', 'Leiden-local'],
	},
	'clean-up': {
		'result_dir': 'clean-up',
		'plot_dir': 'clean_up/',
		'ego_part_algos': "best",
		'ego_params': ['edges'],
		'clean_up_set': 'all',
		'graph_sets': [
			'om',
			'overlap',
			'mu',
			'facebook',
		],
		'plot_graph_sets': [
			'om',
			'mu',
			'facebook',
			'facebook_bar',
		],
		'remove_algo_parts': ['Ego', ' | EdgesScore', 'Infomap + Surprise | ', 'Leiden + Infomap | ',
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {},
		'plots': ['metrics', 'comm_sizes', 'timings', 'num_comms'],
		'plot_algo_set': ['Leiden-Info', 'Info-Surprise', 'Leiden-Info-noOSLOM',
		                  'Info-Surprise-noOSLOM'],
	},
	'test': {
		'result_dir': 'test',
		'plot_dir': 'test/',
		'ego_part_algos': 'global',
		'ego_params': ['test'],
		'score_per_egonet': False,
		'graph_sets': [
			'test',
			# 'om',
			# 'mu',
			# 'facebook'
		],
		'remove_algo_parts': ['Ego', ' | No Clean Up', ' | EdgesScore',
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {},
		'plot_graph_sets': [
			'facebook_bar',
		],
		# 'plots': ['ego_net_extend', 'ego_net_x_extend', 'timings'],
		'plots': ['metrics', 'ego_net_cluster', 'comm_sizes'],
		'plot_algo_set': ['all'],
	},
}