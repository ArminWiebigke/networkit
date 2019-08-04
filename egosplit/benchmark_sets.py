
def get_benchmark_configs():
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
	'local-cluster',
]


benchmark_configs = {
	'edges-score': {
		'result_dir': 'edges-score',
		'plot_dir': 'extend/edges/score_norm/',
		'hue': 'Algorithm',
		'ego_part_algos': 'fast',
		'ego_params': ['no-extend', 'edges-score'],
		'score_per_egonet': False,
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
		'hue': 'Algorithm',
		'ego_part_algos': 'fast',
		'ego_params': ['no-extend', 'edges-factor'],
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
		'hue': 'Algorithm',
		'ego_part_algos': 'leiden local',
		'ego_params': ['no-extend', 'sig-merge'],
		'score_per_egonet': False,
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
		'hue': 'Algorithm',
		'ego_part_algos': 'leiden local',
		'ego_params': ['no-extend', 'sig-ext-iter'],
		'score_per_egonet': False,
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
		'hue': 'Algorithm',
		'ego_part_algos': 'leiden local',
		'ego_params': ['no-extend', 'sig-check-updated'],
		'score_per_egonet': False,
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
		'hue': 'Algorithm',
		'ego_part_algos': 'leiden local',
		'ego_params': ['no-extend', 'sig-max-candidates'],
		'score_per_egonet': False,
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
		'hue': 'Algorithm',
		'ego_part_algos': 'leiden local',
		'ego_params': ['no-extend', 'sig-cluster-iter'],
		'score_per_egonet': False,
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
		'hue': 'Algorithm',
		'ego_part_algos': 'leiden local',
		'ego_params': ['no-extend', 'extend'],
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
		'hue': 'Algorithm',
		'ego_part_algos': 'local',
		'ego_params': ['no-extend', 'extend'],
		'score_per_egonet': True,
		'graph_sets': [
			'om',
			'overlap',
			'mu',
			'facebook',
			'facebook_bar'
		],
		'remove_algo_parts': ['Ego', ' + Infomap', ' | No Clean Up', '!', ' | No Extension',
		                      ' | Edges', ' | Significance'
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {},
		'plots': ['ego_net_cluster', 'metrics', 'timings'],
		'plot_algo_set': ['base', 'edges', 'sig'],
	},
	# 'ext-compare': {
	# 	'result_dir': '*15:18:42',
	# 	'plot_dir': 'extend/comparison/',
	# 	'hue': 'Algorithm',
	# 	'graph_sets': ['om', 'mu', 'facebook'],
	# 	'remove_algo_parts': ['Ego', 'Leiden', ' + Infomap_', '_no-clean', '!',
	# 	                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)]
	# },
	# 'sig-max-checked-cand': {
	# 	'result_dir': '*-sig-max-checked-cand',
	# 	'plot_dir': "extend/sig/max_candidates_checked/",
	# 	'hue': 'Algorithm',
	# 	'graph_sets': ['om'],
	# 	'remove_algo_parts': ['Ego', 'Leiden', ' + Infomap_', '_no-clean', '!',
	# 	                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)]
	# },
	# 'glob-compare-ext': {
	# 	'result_dir': '*-glob-compare-ext',
	# 	'plot_dir': "global_compare_ext/",
	# 	'hue': 'Algorithm',
	# 	'graph_sets': ['om'],
	# 	'remove_algo_parts': ['Ego', '_no-clean', '!',
	# 	                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)]
	# },
	# 'clean-up': {
	# 	'result_dir': '*-clean-up*',
	# 	'plot_dir': "clean_up/",
	# 	'hue': 'Algorithm',
	# 	'graph_sets': ['om', 'mu', 'facebook'],
	# 	'remove_algo_parts': ['Ego', 'map' '!',
	# 	                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)]
	# },
	'test': {
		'result_dir': 'test',
		'plot_dir': 'test/',
		'hue': 'Algorithm',
		'ego_part_algos': 'leiden local',
		'ego_params': ['no-extend', 'test'],
		'score_per_egonet': True,
		'graph_sets': [
			'om',
			# 'mu',
			# 'facebook'
		],
		'remove_algo_parts': ['Ego', 'Leiden + Infomap | ', ' | No Clean Up', '!',
		                      ] + ['({:03.0f})'.format(i) for i in range(0, 30)],
		'replace_legend': {},
		# 'plots': ['ego_net_extend', 'ego_net_x_extend', 'timings'],
		'plots': ['metrics'],
		'plot_algo_set': ['all'],
	},
}
