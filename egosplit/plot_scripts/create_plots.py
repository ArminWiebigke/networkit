from matplotlib.ticker import MultipleLocator

from .config import metric_names
from .draw_plot import make_plot, PlotType
import egosplit.benchmarks.evaluation.benchmark_metric as bm

# algo_matches = [
# '',
# 'Ego',
# 	# '_b#',
# 	'_e#',
# 	# 'e-score',
# 	# 'Leiden',
# 	# 'GCE',
# 	# 'MOSES',
# 	# 'OSLOM'
# ]

# TODO: Extract algorithm parameter as column
remove_algo_part = [
	'Ego_',
	# 'PLM_',
	'Leiden_Mod_',
	'ean-merge,overl',
	'ean',
	'_remv-overl',
	'!',
	'_no-cl',
	'groups',
	# '_b'
]

# TODO: metric: NMI for partition vs. ground-truth, ignore all external nodes
# TODO: metric: for each community, check if the largest part dominates its partition.
#
ego_metrics = [
	'Community Segmentation',
	'Community Merging',
	# 'Merged External Nodes',
	# 'Cluster per Community',
	# 'Communities per Cluster',
	'Community Coverage',
	'Extended Nodes / Ego-Net Size',
	'External Nodes Ratio',
	'Added External Nodes / Added Nodes',
	'Persona Recall',
	# 'Conductance',
	'Community Fitness (alpha=0.8)',
	'Community Fitness (alpha=0.9)',
	'Community Fitness (alpha=1.0)',
	'Community Fitness (alpha=1.1)',
	# 'conductance_ratio',
	# 'Ratio of Intra-Community Edges',
	'Change of Ratio of Intra-Community Edges',
	'Components per Community',
	'Ratio of Disconnected Nodes',
]


def safe_filename(filename):
	return filename.replace("/", "_")


algo_sets = {
	# 'base': ['b#'],
	# 'edges': ['e#'],
	# 'signif': ['s#'],
	'all': ['Ego'],
}

graph_sets = {
	'om': {
		'graph_filter': '_om_',
		'x_filter': None,
		'x': 'Communities per Node',
		'ax_set': {
		}
	},
	# 'all': {
	# 	'graph_filter': '',
	# 	'x_filter': None,
	# 	'x': 'Graph Name',
	# 	'ax_set': {
	# 	},
	# 	'bar_plot': True,
	# },
	# 'mu': {
	# 	'graph_filter': '_mu_',
	# 	'x': 'Mixing Factor',
	# 	'ax_set': {
	# 	}
	# },
	# 'Facebook': {
	# 	'graph_filter': 'FB_',
	# 	'x': 'Graph Name',
	# 	'ax_set': {
	# 	},
	# 	'bar_plot': True,
	# },
}

timings = [
	'___Partition_EgoNet',
	'___Extend_EgoNet',
	'___Extend_and_Partition_EgoNet',
]


def run(data, output_dir):
	plot_funcs = [
		metric_plots,
		timings_plots,
		num_comms_plots,
		ego_net_plots,
		# ego_net_plots_per_graph,
		# comm_sizes_plots,
	]
	for plot_func in plot_funcs:
		for graph_set_name, graph_set_params in graph_sets.items():
			for algo_set_name, algo_match in algo_sets.items():
				params = (data, output_dir, graph_set_name, graph_set_params, algo_set_name,
				          algo_match)
				plot_func(*params)


# *****************************************************************************
# *                                  Timings                                  *
# *****************************************************************************
def timings_plots(data, output_dir, graph_set_name, graph_set_params, algo_set_name, algo_match):
	print('Plots for Timings')
	for timing in timings:
		make_plot(
			output_dir=output_dir,
			plot_type=(PlotType.bar if 'bar_plot' in graph_set_params else PlotType.line),
			data=data['timings'],
			graph_filter=graph_set_params['graph_filter'],
			x_filter=graph_set_params['x_filter'],
			filter_data='`Timer Name`.str.contains(\'{}\')'.format(timing),
			algo_matches=algo_match,
			remove_algo_part=remove_algo_part,
			# title=timing.strip('_'),
			file_name='timings/{}_{}_{}'.format(timing.strip('_'), graph_set_name, algo_set_name),
			one_plot_per_graph=False,
			x=graph_set_params['x'],
			y='Running Time / Average Degree',
			hue='Algorithm',
			plot_args={
				'ci': 'sd',
			},
			ax_set={
				**graph_set_params['ax_set'],
				'ylim': 0,
				'ylabel': r'Running Time / k [ms]',
			}
		)


# *****************************************************************************
# *                                  Metrics                                  *
# *****************************************************************************
def metric_plots(data, output_dir, graph_set_name, graph_set_params, algo_set_name, algo_match):
	print('Plots for Metrics')
	metrics = [
		bm.NMI.get_name(),
		bm.Time.get_name(),
		bm.F1.get_name(),
		bm.F1_rev.get_name(),
	]
	for metric in metrics:
		make_plot(
			output_dir=output_dir,
			plot_type=(PlotType.bar if 'bar_plot' in graph_set_params else PlotType.line),
			data=data['metrics'],
			graph_filter=graph_set_params['graph_filter'],
			x_filter=graph_set_params['x_filter'],
			algo_matches=algo_match,
			remove_algo_part=remove_algo_part,
			# title=metric,
			file_name='metrics/{}_{}_{}'.format(metric,
			                                    graph_set_name, algo_set_name),
			one_plot_per_graph=False,
			x=graph_set_params['x'],
			y=metric,
			hue='Algorithm',
			plot_args={
				'ci': 'sd',
			},
			ax_set={
				**graph_set_params['ax_set'],
				'ylim': metric_names[metric]['ylim'],
				'ylabel': metric,
			}
		)


# *****************************************************************************
# *                                  Comm sizes                               *
# *****************************************************************************
def comm_sizes_plots(data, output_dir, graph_set_name, graph_set_params, algo_set_name, algo_match):
	print('Plots for comm sizes')
	make_plot(
		output_dir=output_dir,
		plot_type=PlotType.swarm,
		data=data['cover_comm_sizes'],
		graph_filter=graph_set_params['graph_filter'],
		algo_matches=algo_match,
		add_algos=['Ground_Truth'],
		remove_algo_part=remove_algo_part,
		# title='Community Sizes',
		file_name='communities/comm_sizes_{}_{}'.format(graph_set_name, algo_set_name),
		one_plot_per_graph=True,
		x='Graph Name',
		hue='Algorithm',
		y='Community Size',
		plot_args={
			# 'size': 2.5,
			'dodge': True,
		},
		ax_set={
			# 'ylim': (2, 8),
			# 'ylim': 2,
			'ylabel': 'size (log2)',
			# 'xticklabels': [''],
		}
	)


# *****************************************************************************
# *                            Number of communities                          *
# *****************************************************************************
def num_comms_plots(data, output_dir, graph_set_name, graph_set_params, algo_set_name, algo_match):
	print('Plots for num comms')
	make_plot(
		output_dir=output_dir,
		# plot_type=PlotType.bar,
		data=data['cover_num_comms'],
		graph_filter=graph_set_params['graph_filter'],
		x_filter=graph_set_params['x_filter'],
		algo_matches=algo_match,
		add_algos=['Ground_Truth'],
		remove_algo_part=remove_algo_part,
		# title='Number of communities',
		file_name='communities/num_comms_{}_{}'.format(graph_set_name, algo_set_name),
		one_plot_per_graph=False,
		x=graph_set_params['x'],
		y='Number of Communities',
		hue='Algorithm',
		plot_args={
		},
		ax_set={
			'ylim': 0,
		}
	)


# *****************************************************************************
# *                             Ego-Net partition                             *
# *****************************************************************************
def ego_net_plots(data, output_dir, graph_set_name, graph_set_params, algo_set_name, algo_match):
	print('Plots for ego-net metrics')
	for ego_metric in ego_metrics:
		make_plot(
			output_dir=output_dir,
			plot_type=(PlotType.bar if 'bar_plot' in graph_set_params else PlotType.line),
			data=data['ego_net_metrics'],
			filter_data='Metric == \'{}\''.format(ego_metric),
			graph_filter=graph_set_params['graph_filter'],
			x_filter=graph_set_params['x_filter'],
			algo_matches=algo_match,
			remove_algo_part=remove_algo_part,
			# title='Ego-Net {}'.format(ego_metric),
			file_name='ego_partition/metrics/{}_{}_{}'.format(safe_filename(ego_metric),
			                                                  graph_set_name,
			                                                  algo_set_name),
			one_plot_per_graph=False,
			x=graph_set_params['x'],
			y='Value',
			hue='Algorithm',
			plot_args={
				# 'style': 'metric_name',
				'ci': 'sd',
			},
			ax_set={
				# 'ylim': (0, 1.05),
				'ylim': 0,
				'ylabel': ego_metric,
			}
		)


# # Algo parameters on x-axis
# for ego_metric in ego_metrics:
# 	break
# 	x = {'name': 'factor', 'create_from': 'Algorithm', 'str_start': '_f-', 'str_end': '*'}
# 	hue = {'name': 'exponent', 'create_from': 'Algorithm', 'str_start': '*e-',
# 	       'str_end': ''}
# 	make_plot(
# 		output_dir=output_dir,
# 		data=data['ego_net_metrics'].query('metric_name in @ego_metric'),
# 		graphs='LFR_om',
# 		algo_match=algo_matches,
# 		title='Ego-Net Metrics, ' + ego_metric,
# 		file_name='ego_partition/2_dim_metrics/' + ego_metric,
# 		x=x,
# 		y='value',
# 		hue=hue,
# 		plot_args={
# 		},
# 		ax_set={
# 			'ylim': (0, 1.05),
# 		}
# 	)


def ego_net_plots_per_graph(data, output_dir, graph_set_name, graph_set_params, algo_set_name,
                            algo_match):
	# Metrics per Ego-Net
	print('Plots for ego-net metrics per graph')
	for ego_metric in ego_metrics:
		make_plot(
			output_dir=output_dir,
			# plot_type=PlotType.bar,
			data=data['ego_net_ego_metrics'],
			filter_data=('Metric in \'{}\''.format(ego_metric)),
			graph_filter=graph_set_params['graph_filter'],
			algo_matches=algo_match,
			remove_algo_part=remove_algo_part,
			# title='Ego-Net {}'.format(ego_metric),
			file_name='ego_partition/ego_metrics/{}_{}_{}'.format(safe_filename(ego_metric),
			                                                      graph_set_name,
			                                                      algo_set_name),
			one_plot_per_graph=True,
			x='Ego-Net Size',
			y='Value',
			hue='Algorithm',
			plot_args={
				# 'style': 'metric_name',
				# 'markersize': 3,
			},
			ax_set={
				# 'ylim': (0, 1.05),
				'ylim': 0,
				'ylabel': ego_metric,
			}
		)
