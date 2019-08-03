from collections import defaultdict

from matplotlib.ticker import MultipleLocator

from .config import metric_names
from .draw_plot import make_plot, PlotType
import egosplit.benchmarks.evaluation.benchmark_metric as bm

ego_extend_metrics = {
	'Community Coverage': 'coverage',
	'Extended Nodes / Ego-Net Size': 'extended_nodes',
	'External Nodes Ratio': "external_nodes",
	'Added External Nodes Ratio': 'added_external',
	# 'Conductance': 'conductance',
	'Community Fitness': 'comm_fitness',
	'Ratio of Intra-Community Edges': 'intra_ratio',
	# 'Change of Intra-Edges': 'intra_change',
	'Components per Community': 'components',
	'Ratio of Disconnected Nodes': 'disconnected_nodes',
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

algo_sets = {
	# 'base': ['b!'],
	# 'edges': ['e!'],
	# 'signif': ['s!'],
	'all': ['Ego', 'GCE', 'OSLOM', 'MOSES'],
	# 'ego': ['Ego'],
}

graph_sets = {
	'om': {
		'graph_filter': '_om_',
		'x': 'Communities per Node',
		'x_filter': None,
		'ax_set': {
		}
	},
	'all': {
		'graph_filter': '',
		'x': 'Graph Name',
		'x_filter': None,
		'ax_set': {
		},
		'bar_plot': True,
	},
	'mu': {
		'graph_filter': '_mu_',
		'x': 'Mixing Factor',
		'x_filter': None,
		'ax_set': {
		}
	},
	'facebook': {
		'graph_filter': 'FB_',
		'x': 'Graph Name',
		'x_filter': None,
		'plot_args': {
			'dashes': [(1, 2) for _ in range(10)],
		},
		'ax_set': {
		},
		# 'bar_plot': True,
	},
}

timings = [
	'   Partition EgoNet',
	'   Extend EgoNet',
	'   Extend and Partition EgoNet',
]

metrics = [
	bm.NMI.get_name(),
	bm.Time.get_name(),
	bm.F1.get_name(),
	bm.F1_rev.get_name(),
]


def run(data, output_dir, config):
	assert output_dir[-1] == "/"
	plot_funcs = {
		'metrics': metric_plots,
		'timings': timings_plots,
		'ego_net_extend': ego_net_extend_plots,
		'ego_net_cluster': ego_net_cluster_plots,
		'ego_net_x_extend': ego_net_x_extend_plots,
		'ego_net_x_cluster': ego_net_x_cluster_plots,
		'num_comms': num_comms_plots,
		'comm_sizes': comm_sizes_plots,
	}
	for plot_func in config['plots']:
		for graph_set_name in config['graph_sets']:
			graph_set_params = graph_sets[graph_set_name]
			algo_set_name = config['plot_algo_set']
			algo_match = algo_sets[algo_set_name]
			params = (data, output_dir, graph_set_name, graph_set_params, algo_set_name,
			          algo_match, config)
			plot_funcs[plot_func](*params)


def safe_filename(filename):
	return filename.replace("/", "_")


# *****************************************************************************
# *                                  Timings                                  *
# *****************************************************************************
def timings_plots(data, output_dir, graph_set_name, graph_set_params, algo_set_name, algo_match,
                  config):
	print('Plots for Timings')
	for timing in timings:
		make_plot(
			plot_subdir='timings/',
			output_dir=output_dir,
			plot_type=(PlotType.bar if 'bar_plot' in graph_set_params else PlotType.line),
			data=data['timings'],
			graph_filter=graph_set_params['graph_filter'],
			x_filter=graph_set_params['x_filter'],
			filter_data='`Timer Name`.str.contains(\'{}\')'.format(timing),
			algo_matches=algo_match,
			remove_algo_part=config['remove_algo_parts'],
			replace_legend=config['replace_legend'],
			# title=timing.strip(),
			file_name='{}_{}_{}'.format(timing.strip(' '), graph_set_name, algo_set_name),
			legend_file_name='legend_{}_{}'.format(graph_set_name, algo_set_name),
			one_plot_per_graph=False,
			x=graph_set_params['x'],
			y='Running Time / (n + m)',
			# y='Running Time / Average Degree',
			hue=config['hue'],
			plot_args={
				**graph_set_params.get('plot_args', {}),
				# **{k: v for k, v in graph_set_params.items() if k == 'line_dashes'},
				'ci': 'sd' if config.get('show_deviation', False) else None,
			},
			ax_set={
				**graph_set_params['ax_set'],
				'ylim': 0,
				'ylabel': r'Running Time / (n + m) [\SI{}{\micro\second}]',
			}
		)


# *****************************************************************************
# *                                  Metrics                                  *
# *****************************************************************************
def metric_plots(data, output_dir, graph_set_name, graph_set_params, algo_set_name, algo_match,
                 config):
	print('Plots for Metrics')
	for metric in metrics:
		make_plot(
			plot_subdir='metrics/',
			output_dir=output_dir,
			plot_type=(PlotType.bar if 'bar_plot' in graph_set_params else PlotType.line),
			data=data['metrics'],
			graph_filter=graph_set_params['graph_filter'],
			x_filter=graph_set_params['x_filter'],
			algo_matches=algo_match,
			remove_algo_part=config['remove_algo_parts'],
			# title=metric,
			file_name='{}_{}_{}'.format(metric, graph_set_name, algo_set_name),
			legend_file_name='legend_{}_{}'.format(graph_set_name, algo_set_name),
			one_plot_per_graph=False,
			x=graph_set_params['x'],
			y=metric,
			hue=config['hue'],
			plot_args={
				'ci': 'sd' if config.get('show_deviation', False) else None,
			},
			ax_set={
				**graph_set_params['ax_set'],
				'ylim': metric_names[metric]['ylim'],
				'ylabel': metric_names[metric]['ylabel'],
			}
		)


# *****************************************************************************
# *                                  Comm sizes                               *
# *****************************************************************************
def comm_sizes_plots(data, output_dir, graph_set_name, graph_set_params, algo_set_name, algo_match,
                     config):
	print('Plots for comm sizes')
	make_plot(
		plot_subdir='comm_sizes/',
		output_dir=output_dir,
		plot_type=PlotType.swarm,
		data=data['cover_comm_sizes'],
		graph_filter=graph_set_params['graph_filter'],
		algo_matches=algo_match,
		add_algos=['Ground Truth'],
		remove_algo_part=config['remove_algo_parts'],
		# title='Community Sizes',
		file_name='comm_sizes_{}_{}'.format(graph_set_name, algo_set_name),
		legend_file_name='legend_{}_{}'.format(graph_set_name, algo_set_name),
		one_plot_per_graph=True,
		x='Graph Name',
		hue=config['hue'],
		y='Community Size',
		plot_args={
			# 'size': 2.5,
			'dodge': True,
		},
		ax_set={
			# 'ylim': (2, 8),
			# 'ylim': 2,
			'ylabel': 'Community Size ($\log_2$)',
			# 'xticklabels': [''],
		},
		# font_size=0.6,
	)


# *****************************************************************************
# *                            Number of communities                          *
# *****************************************************************************
def num_comms_plots(data, output_dir, graph_set_name, graph_set_params, algo_set_name, algo_match,
                    config):
	print('Plots for num comms')
	make_plot(
		plot_subdir='num_comms/',
		output_dir=output_dir,
		# plot_type=PlotType.bar,
		data=data['cover_num_comms'],
		graph_filter=graph_set_params['graph_filter'],
		x_filter=graph_set_params['x_filter'],
		algo_matches=algo_match,
		add_algos=['Ground Truth'],
		remove_algo_part=config['remove_algo_parts'],
		# title='Number of communities',
		file_name='num_comms_{}_{}'.format(graph_set_name, algo_set_name),
		legend_file_name='legend_{}_{}'.format(graph_set_name, algo_set_name),
		one_plot_per_graph=False,
		x=graph_set_params['x'],
		y='Number of Communities',
		hue=config['hue'],
		plot_args={
		},
		ax_set={
			'ylim': 0,
		}
	)


# *****************************************************************************
# *                             Ego-Net partition                             *
# *****************************************************************************
def ego_net_extend_plots(*args):
	ego_net_plots(*args, ego_extend_metrics)


def ego_net_cluster_plots(*args):
	ego_net_plots(*args, ego_cluster_metrics)


def ego_net_plots(data, output_dir, graph_set_name, graph_set_params, algo_set_name, algo_match,
                  config, ego_metrics):
	print('Plots for ego-net metrics')
	for ego_metric, file_name in ego_metrics.items():
		make_plot(
			plot_subdir='ego_metrics/',
			output_dir=output_dir,
			plot_type=(PlotType.bar if 'bar_plot' in graph_set_params else PlotType.line),
			data=data['ego_net_metrics'],
			filter_data='Metric == \'{}\''.format(ego_metric),
			graph_filter=graph_set_params['graph_filter'],
			x_filter=graph_set_params['x_filter'],
			algo_matches=algo_match,
			remove_algo_part=config['remove_algo_parts'],
			replace_legend=config['replace_legend'],
			# title='Ego-Net {}'.format(ego_metric),
			file_name='{}_{}_{}'.format(file_name, graph_set_name, algo_set_name),
			legend_file_name='legend_{}_{}'.format(graph_set_name, algo_set_name),
			one_plot_per_graph=False,
			x=graph_set_params['x'],
			y='Value',
			hue=config['hue'],
			plot_args={
				# 'style': 'metric_name',
				'ci': 'sd' if config.get('show_deviation', False) else None,
			},
			ax_set={
				# 'ylim': (0, 1.05),
				'ylim': 0,
				'ylabel': ego_metric,
			}
		)


def ego_net_x_extend_plots(*args):
	ego_net_plots_per_graph(*args, ego_extend_metrics)


def ego_net_x_cluster_plots(*args):
	ego_net_plots_per_graph(*args, ego_cluster_metrics)


def ego_net_plots_per_graph(data, output_dir, graph_set_name, graph_set_params, algo_set_name,
                            algo_match, config, ego_metrics):
	# Metrics per Ego-Net
	print('Plots for ego-net metrics per graph')
	for ego_metric, file_name in ego_metrics.items():
		make_plot(
			plot_subdir='ego_size_metrics/',
			output_dir=output_dir,
			# plot_type=PlotType.bar,
			data=data['ego_net_ego_metrics'],
			filter_data=('Metric in \'{}\''.format(ego_metric)),
			graph_filter=graph_set_params['graph_filter'],
			algo_matches=algo_match,
			remove_algo_part=config['remove_algo_parts'],
			# title='Ego-Net {}'.format(ego_metric),
			file_name='{}_{}_{}'.format(file_name, graph_set_name, algo_set_name),
			legend_file_name='legend_{}_{}'.format(graph_set_name, algo_set_name),
			one_plot_per_graph=True,
			x='Ego-Net Size',
			y='Value',
			hue=config['hue'],
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
