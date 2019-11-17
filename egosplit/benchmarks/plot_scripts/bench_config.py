from collections import defaultdict

from egosplit.benchmarks.evaluation.metrics import NMI, F1, F1_rev, Time

ego_extend_metrics = {
	'Community Coverage': 'coverage',
	'Extended Nodes / Ego-Net Size': 'extended_nodes',
	'External Nodes Ratio': 'external_nodes',
	'Added External Nodes Ratio': 'added_external',
	# 'Conductance': 'conductance',
	'Community Fitness': 'comm_fitness',
	'Ratio of Intra-Community Edges': 'intra_ratio',
	# 'Change of Intra-Edges': 'intra_change',
	'Components per Community': 'components',
	'Ratio of Disconnected Nodes': 'disconnected_nodes',
	'Number of Ground-Truth Communities': 'gt_comms',
}

ego_cluster_metrics = {
	'Community Segmentation': 'segmentation',
	'Community Merging': 'comm_merge',
	'Merged External Nodes': 'merged_external',
	'Cluster per Community': 'cluster_per_comm',
	'Communities per Cluster': 'comms_per_cluster',
	'Number of Personas': 'personas',
	'Number of Ground-Truth Communities': 'gt_comms',
	'Persona Recall': 'persona_recall',
}


class PlotAlgoSetConfig:
	@staticmethod
	def get_algo_sets(algo_set_names):
		algo_set_dict = {
			'base': ['No Extension'],
			'edges': ['Edges'],
			'sig': ['Significance'],
			'all': ['!Ground Truth', ''],
			'ego': ['Ego'],
			'connect-persona': ['Leiden + Infomap', '!All Unweighted'],
			'clean-up': ['Leiden + Infomap', '!Clean-remove', '!OSLOM-full'],
			'Info-Info': ['Infomap + Infomap'],
			'Info-local': ['Infomap + '],
			'Leiden-local': ['Leiden + '],
			'MapEquation-local': ['MapEquation + '],
			'Info-Surprise': ['Infomap + Surprise', '!Remove Overlapping'],
			'Info-Surprise-noOSLOM': ['Infomap + Surprise', '!OSLOM-full', '!Remove Overlapping'],
			'Leiden-Info': ['Leiden + Infomap', '!Remove Overlapping'],
			'Leiden-Info-noOSLOM': ['Leiden + Infomap', '!OSLOM-full', '!Remove Overlapping'],
			'PLP': ['PLP'],
			'PLM': ['PLM'],
			'Potts': ['Potts'],
		}
		algo_sets = []
		for algo_set_name in algo_set_names:
			algo_sets.append((algo_set_name, algo_set_dict[algo_set_name]))
		return algo_sets


algo_sets = {
	'base': ['No Extension'],
	'edges': ['Edges'],
	'sig': ['Significance'],
	'all': ['!Ground Truth', ''],
	'ego': ['Ego'],
	'connect-persona': ['Leiden + Infomap', '!All Unweighted'],
	'clean-up': ['Leiden + Infomap', '!Clean-remove', '!OSLOM-full'],
	'Info-Info': ['Infomap + Infomap'],
	'Info-local': ['Infomap + '],
	'Leiden-local': ['Leiden + '],
	'Info-Surprise': ['Infomap + Surprise', '!Remove Overlapping'],
	'Info-Surprise-noOSLOM': ['Infomap + Surprise', '!OSLOM-full', '!Remove Overlapping'],
	'Leiden-Info': ['Leiden + Infomap', '!Remove Overlapping'],
	'Leiden-Info-noOSLOM': ['Leiden + Infomap', '!OSLOM-full', '!Remove Overlapping'],
}
for key in algo_sets:
	algo_sets[key] = (key, algo_sets[key])

graph_sets = {
	'all': {
		'name': 'all',
		'graph_filter': '',
		'x': 'Graph Name',
		'x_filter': None,
		'ax_set': {
		},
		'bar_plot': True,
	},
	'om': {
		'name': 'om',
		'graph_filter': '_om_',
		'x': 'Communities per Node',
		'x_filter': None,
		'ax_set': {
		}
	},
	'om_max': {
		'name': 'om_max',
		'graph_filter': '_om_',
		'x': 'Communities per Node',
		'x_filter': None,
		'ax_set': {
			'ylim': (0, 1190),
		}
	},
	'mu': {
		'name': 'mu',
		'graph_filter': '_mu_',
		'x': 'Mixing Factor',
		'x_filter': None,
		'ax_set': {
		}
	},
	'mu_max': {
		'name': 'mu_max',
		'graph_filter': '_mu_',
		'x': 'Mixing Factor',
		'x_filter': None,
		'ax_set': {
			'ylim': (0, 990),
		}
	},
	'om_bar': {
		'name': 'om_bar',
		'graph_filter': '_om_',
		'x': 'Communities per Node',
		'x_filter': None,
		'ax_set': {
		},
		'bar_plot': True,
	},
	'facebook': {
		'name': 'facebook',
		'graph_filter': 'FB_',
		'x': 'Graph Name',
		'x_filter': None,
		'plot_args': {
			'dashes': [(5, 6) for _ in range(10)],
		},
		'ax_set': {
		},
		'set_ylim': False,
	},
	'facebook_bar': {
		'name': 'facebook_bar',
		'graph_filter': 'FB_',
		'x': 'Graph Name',
		'x_filter': None,
		'plot_args': {
		},
		'ax_set': {
		},
		'bar_plot': True,
		'set_ylim': False,
		'show_deviation': True,
	},
	'large': {
		'name': 'large',
		'graph_filter': ['Amazon', 'DBLP'],
		'x': 'Graph Name',
		'x_filter': None,
		'plot_args': {
		},
		'ax_set': {
		},
		'bar_plot': True,
		'set_ylim': False,
	},
	'overlap': None
}


class PlotGraphSetConfig:
	@staticmethod
	def get_sets(graph_set_names):
		graph_set_dict = {
			'all': {
				'name': 'all',
				'graph_filter': '',
				'x': 'Graph Name',
				'x_filter': None,
				'ax_set': {
				},
				'bar_plot': True,
			},
			'om': {
				'name': 'om',
				'graph_filter': '_om_',
				'x': 'Communities per Node',
				'x_filter': None,
				'ax_set': {
				}
			},
			'om_max': {
				'name': 'om_max',
				'graph_filter': '_om_',
				'x': 'Communities per Node',
				'x_filter': None,
				'ax_set': {
					'ylim': (0, 1190),
				}
			},
			'mu': {
				'name': 'mu',
				'graph_filter': '_mu_',
				'x': 'Mixing Factor',
				'x_filter': None,
				'ax_set': {
				}
			},
			'mu_max': {
				'name': 'mu_max',
				'graph_filter': '_mu_',
				'x': 'Mixing Factor',
				'x_filter': None,
				'ax_set': {
					'ylim': (0, 990),
				}
			},
			'om_bar': {
				'name': 'om_bar',
				'graph_filter': '_om_',
				'x': 'Communities per Node',
				'x_filter': None,
				'ax_set': {
				},
				'bar_plot': True,
			},
			'facebook': {
				'name': 'facebook',
				'graph_filter': 'FB_',
				'x': 'Graph Name',
				'x_filter': None,
				'plot_args': {
					'dashes': [(5, 6) for _ in range(10)],
				},
				'ax_set': {
				},
				'set_ylim': False,
			},
			'facebook_bar': {
				'name': 'facebook_bar',
				'graph_filter': 'FB_',
				'x': 'Graph Name',
				'x_filter': None,
				'plot_args': {
				},
				'ax_set': {
				},
				'bar_plot': True,
				'set_ylim': False,
				'show_deviation': True,
			},
			'snap': {
				'name': 'snap',
				'graph_filter': ['Amazon', 'DBLP', 'LiveJournal', 'Orkut'],
				'x': 'Graph Name',
				'x_filter': None,
				'plot_args': {
				},
				'ax_set': {
				},
				'bar_plot': True,
				'set_ylim': False,
			},
			'overlap': None
		}
		graph_sets = []
		for graph_set_name in graph_set_names:
			graph_sets.append(graph_set_dict[graph_set_name])
		return graph_sets

timings = [
	'   Partition EgoNet',
	'   Extend EgoNet',
	'   Extend and Partition EgoNet',
	'  create EgoNets',
	'  Persona Clustering',
	' Connect Personas',
]

metrics = [
	NMI.get_name(),
	Time.get_name(),
	F1.get_name(),
	F1_rev.get_name(),
]

metric_config = {
	F1.get_name(): {
		'ylabel': 'F1 Score',
		'y_val': 'F1-Score',
		'file_name': 'F1',
		'ylim': (0, 1.05),
	},
	F1_rev.get_name(): {
		'ylabel': 'F1 Score (reversed)',
		'y_val': 'F1-Score (reversed)',
		'file_name': 'F1_rev',
		'ylim': (0, 1.05),
	},
	NMI.get_name(): {
		'ylabel': 'NMI',
		'y_val': 'NMI',
		'file_name': 'NMI',
		'ylim': (0, 1.05),
	},
	Time.get_name(): {
		'ylabel': r'Running Time / (n + m) [\SI{}{\micro\second}]',
		'y_val': 'time / (n + m)',
		'file_name': 'time',
		'ylim': 0,
	},
}

ego_metric_ylim = defaultdict(lambda: 0)
ego_metric_ylim['Persona Recall'] = (0, 1.05)
