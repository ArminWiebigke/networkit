from egosplit.benchmarks.plot_scripts.make_plot import make_plot, PlotType
from egosplit.benchmarks.plot_scripts.bench_config import ego_extend_metrics, ego_cluster_metrics, \
	algo_sets, timings, metrics, metric_config, ego_metric_ylim, PlotGraphSetConfig, \
	PlotAlgoSetConfig


class PlotSetConfig:
	@staticmethod
	def get_plot_functions(plot_sets):
		function_dict = {
			'metrics': metric_plots,
			'timings': timings_plots,
			'comm_sizes': comm_sizes_plots,
			'comm_f1': comm_f1_plots,
			'num_comms': num_comms_plots,
			'ego_net_extend': ego_net_extend_plots,
			'ego_net_cluster': ego_net_cluster_plots,
			'ego_net_x_extend': ego_net_x_extend_plots,
			'ego_net_x_cluster': ego_net_x_cluster_plots,
		}
		plot_functions = []
		for plot_set in plot_sets:
			plot_functions.append(function_dict[plot_set])
		return plot_functions


def make_plots(config, data, algo_set_name, algo_set, graph_set_params, output_dir, plot_func):
	print(algo_set)
	print(graph_set_params)
	graph_set_name = graph_set_params['name']
	default_plots_config = {
		'output_dir': output_dir,
		'plot_type': (
			PlotType.bar if 'bar_plot' in graph_set_params else PlotType.line),
		'graph_filter': graph_set_params['graph_filter'],
		'x_filter': graph_set_params['x_filter'],
		'algo_matches': algo_set,
		'remove_algo_part': config['remove_algo_parts'],
		'replace_legend': config['replace_legend'],
		'legend_file_name': 'legend_{}_{}'.format(graph_set_name, algo_set_name),
		'x': graph_set_params['x'],
		'hue': config['hue'],
		'ax_set': graph_set_params.get('ax_set', {}),
		'plot_args': graph_set_params.get('plot_args', {}),
	}
	plot_func(
		default_plots_config,
		data=data,
		config=config,
		graph_set_name=graph_set_name,
		graph_set_params=graph_set_params,
		algo_set_name=algo_set_name,
		output_dir=output_dir,
	)


def metric_plots(plots_config, **kwargs):
	print('Plots for Metrics')
	for metric in metrics:
		metric_plots_config = {
			'plot_subdir': 'metrics/',
			'data': kwargs['data']['metrics'],
			'file_name': '{}_{}_{}'.format(metric, kwargs['graph_set_name'],
			                               kwargs['algo_set_name']),
			'y': metric_config[metric]['y_val'],
			'ax_set': {
				**plots_config['ax_set'],
				'ylim': metric_config[metric]['ylim'] if kwargs['graph_set_params'].get('set_ylim',
				                                                                        True) else 0,
				'ylabel': metric_config[metric]['ylabel'],
			}
		}
		metric_plots_config = {**plots_config, **metric_plots_config}
		make_plot(**metric_plots_config)


def timings_plots(default_plots_config, **kwargs):
	print('Plots for Timings')
	for timing in timings:
		timings_plots_config = {
			'plot_subdir': 'timings/',
			'data': kwargs['data']['timings'],
			'filter_data': '`Timer Name`.str.contains(\'{}\')'.format(timing),
			'file_name': '{}_{}_{}'.format(timing.strip(' '), kwargs['graph_set_name'],
			                               kwargs['algo_set_name']),
			'y': 'Running Time / (n + m)',
			'ax_set': {
				**default_plots_config['ax_set'],
				'ylim': 0,
				'ylabel': r'Running Time / (n + m) [\SI{}{\micro\second}]',
			}
		}
		timings_plots_config = {**default_plots_config, **timings_plots_config}
		make_plot(**timings_plots_config)


def comm_sizes_plots(default_plots_config, **kwargs):
	print('Plots for comm sizes')
	comm_sizes_plots_config = {
		'plot_subdir': 'comm_sizes/',
		'plot_type': PlotType.swarm,
		'data': kwargs['data']['cover_comm_sizes'],
		'add_algos': ['Ground Truth'],
		'file_name': 'comm_sizes_{}_{}'.format(kwargs['graph_set_name'], kwargs['algo_set_name']),
		'one_plot_per_graph': True,
		'use_graph_id': True,
		'x': 'Graph Name',
		'y': 'Community Size',
		'plot_args': {
			'dodge': True,
		},
		'ax_set': {
			**default_plots_config['ax_set'],
			'ylabel': 'Community Size ($\log_2$)',
		},
	}
	comm_sizes_plots_config = {**default_plots_config, **comm_sizes_plots_config}
	make_plot(**comm_sizes_plots_config)


def comm_f1_plots(default_plots_config, **kwargs):
	print('Plots for comm sizes')

	comm_f1_plots_config = {
		'plot_subdir': 'comm_f1/',
		'plot_type': PlotType.swarm,
		'data': kwargs['data']['cover_comm_sizes'],
		'file_name': '{}_{}'.format(kwargs['graph_set_name'], kwargs['algo_set_name']),
		'one_plot_per_graph': True,
		'use_graph_id': True,
		'x': 'Graph Name',
		'y': 'F1 Score',
		'plot_args': {
			# 'size': 2.5,
			'dodge': True,
		},
		'ax_set': {
			**default_plots_config['ax_set'],
			'ylabel': 'F1 Score',
		}
	}
	comm_f1_plots_config = {**default_plots_config, **comm_f1_plots_config}
	make_plot(**comm_f1_plots_config)


def num_comms_plots(default_plots_config, **kwargs):
	print('Plots for num comms')
	num_comms_plots_config = {
		'plot_subdir': 'num_comms/',
		'data': kwargs['data']['cover_num_comms'],
		'add_algos': ['Ground Truth'],
		'file_name': 'num_comms_{}_{}'.format(kwargs['graph_set_name'], kwargs['algo_set_name']),
		'y': 'Number of Communities',
		'ax_set': {
			**default_plots_config['ax_set'],
			'ylim': 0,
		}
	}
	num_comms_plots_config = {**default_plots_config, **num_comms_plots_config}
	make_plot(**num_comms_plots_config)


def ego_net_extend_plots(default_plots_config, **kwargs):
	ego_net_plots(default_plots_config, ego_extend_metrics, **kwargs)


def ego_net_cluster_plots(default_plots_config, **kwargs):
	ego_net_plots(default_plots_config, ego_cluster_metrics, **kwargs)


def ego_net_plots(default_plots_config, ego_metrics, **kwargs):
	print('Plots for ego-net metrics')
	for ego_metric, file_name in ego_metrics.items():
		ego_net_plots_config = {
			'plot_subdir': 'ego_metrics/',
			'data': kwargs['data']['ego_net_metrics'],
			'filter_data': 'Metric == \'{}\''.format(ego_metric),
			'file_name': '{}_{}_{}'.format(file_name, kwargs['graph_set_name'],
			                               kwargs['algo_set_name']),
			'y': 'Value',
			'ax_set': {
				**default_plots_config['ax_set'],
				# 'ylim': (0, 1.05),
				'ylim': ego_metric_ylim[ego_metric] if kwargs['graph_set_params'].get('set_ylim',
				                                                                      True) else 0,
				'ylabel': ego_metric,
			}
		}
		ego_net_plots_config = {**default_plots_config, **ego_net_plots_config}
		make_plot(**ego_net_plots_config)


def ego_net_x_extend_plots(default_plots_config, **kwargs):
	ego_net_plots_per_graph(default_plots_config, ego_extend_metrics, **kwargs)


def ego_net_x_cluster_plots(default_plots_config, **kwargs):
	ego_net_plots_per_graph(*default_plots_config, ego_cluster_metrics, **kwargs)


def ego_net_plots_per_graph(default_plots_config, ego_metrics, **kwargs):
	# Metrics for different Ego-Net sizes
	print('Plots for ego-net metrics per graph')
	for ego_metric, file_name in ego_metrics.items():
		ego_net_plots_config = {
			'plot_subdir': 'ego_size_metrics/',
			'data': kwargs['data']['ego_net_ego_metrics'],
			'filter_data': ('Metric in \'{}\''.format(ego_metric)),
			'file_name': '{}_{}_{}'.format(file_name, kwargs['graph_set_name'],
			                               kwargs['algo_set_name']),
			'one_plot_per_graph': True,
			'x': 'Ego-Net Size',
			'y': 'Value',
			'plot_args': {
				# 'style': 'metric_name',
				# 'markersize': 3,
			},
			'ax_set': {
				**default_plots_config['ax_set'],
				# 'ylim': (0, 1.05),
				'ylim': ego_metric_ylim[ego_metric] if kwargs['graph_set_params'].get('set_ylim',
				                                                                      True) else 0,
				'ylabel': ego_metric,
			}
		}
		ego_net_plots_config = {**default_plots_config, **ego_net_plots_config}
		make_plot(**ego_net_plots_config)


def safe_filename(filename):
	return filename.replace("/", "_")
