import matplotlib.pyplot as plt
import seaborn as sns

from enum import Enum

from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FormatStrFormatter, \
	StrMethodFormatter, Formatter

from .config import set_layout
from .read_data import create_column_if_missing


class PlotType(Enum):
	line = 1
	bar = 2
	swarm = 3
	violin = 4


# Create a plot
def make_plot(data,
              output_dir,
              filter_data=None,
              x_filter=None,
              graph_filter='',
              algo_matches='',
              add_algos=None,
              remove_algo_part=None,
              one_plot_per_graph=False,
              title='',
              file_name='',
              plot_type=PlotType.line,
              x=None,
              y=None,
              hue=None,
              plot_args=None,
              ax_set=None):
	# Create the columns for x, y, hue if not already present in data
	for column in [x, y, hue]:
		create_column_if_missing(data, column)
	# Filter data
	if isinstance(graph_filter, str):
		filtered_data = data.query('graph.str.contains(@graph_filter)').copy()
	else:
		assert (isinstance(graph_filter, list))
		filtered_data = data.query('graph in @graph_filter').copy()
	if filter_data:
		filtered_data.query(filter_data, inplace=True)
	if x_filter:
		filtered_data.query(x_filter, inplace=True)
	graphs = get_unique_values(filtered_data, 'graph')
	algo_list = get_algo_list(algo_matches, add_algos, filtered_data)
	filtered_data.query('algo in @algo_list', inplace=True)
	if len(filtered_data) is 0:
		return

	# Create plots
	if not output_dir[-1] == "/":
		output_dir += "/"
	def create_plot(graph_data):
		num_x_values = len(get_unique_values(graph_data, x))
		this_plot_type = confirm_plot_type(plot_type, graph_data, num_x_values)
		this_plot_args = get_plot_args(algo_list, hue, plot_args, this_plot_type, x, y)
		fig, ax = plt.subplots()
		this_plot_args = {
			**this_plot_args,
			'data': graph_data,
			'ax': ax,
		}
		draw_plot(this_plot_type, this_plot_args)

		sns.despine(ax=ax)
		set_ax(ax, ax_set, fig, x)
		clean_legend(algo_matches, ax, remove_algo_part)
		fig.suptitle(title)
		fig.savefig(output_dir + file_name + '.pdf')
		plt.close(fig)

	# Make one plot per graph if graph is not on the x-axis
	if one_plot_per_graph:
		file_name_base = file_name
		title_base = title
		for graph in graphs:
			file_name = file_name_base + "_" + graph
			title = "{} ({})".format(title_base, graph)
			graph_data = filtered_data.query('graph == @graph')
			create_plot(graph_data)
	else:
		create_plot(filtered_data)


def confirm_plot_type(plot_type, graph_data, x_values):
	if plot_type == PlotType.swarm and len(graph_data) > 10000:
		return PlotType.violin

	if x_values == 1 and plot_type == PlotType.line:
		plot_type = PlotType.bar
	return plot_type


def get_unique_values(filtered_data, column):
	return filtered_data.groupby(column).mean().index.values


def clean_legend(algo_matches, ax, remove_algo_part):
	legend_handles, legend_labels = ax.get_legend_handles_labels()
	if remove_algo_part:
		assert (isinstance(remove_algo_part, list))
		remove_list = remove_algo_part
	else:
		remove_list = algo_matches
	for remove in remove_list:
		legend_labels = [l.replace(remove, '') for l in legend_labels]
	set_layout(ax, legend_handles, legend_labels)


def draw_plot(plot_type, plot_args):
	plot_functions = {
		PlotType.line: sns.lineplot,
		PlotType.bar: sns.barplot,
		PlotType.swarm: sns.swarmplot,
		PlotType.violin: sns.violinplot,
	}
	plot_functions[plot_type](**plot_args)


class MinorFormatter(Formatter):
	def __call__(self, x, pos=None):
		s = "{:.1f}".format(x)
		return s[-2:]


def set_ax(ax, ax_set, fig, x):
	ax_set = ax_set or {}
	ax.set(
		**ax_set,
	)
	if x == "communities_per_node":
		ax.xaxis.set(
			major_locator=MultipleLocator(1),
			major_formatter=StrMethodFormatter('{x:.0f}'),
			minor_formatter=MinorFormatter(),
		)
		min_x, max_x = ax.get_xlim()
		if min_x <= 1.2 and max_x >= 1.8:
			minor_xticks = [1.2, 1.4, 1.6, 1.8]
			ax.set_xticks(minor_xticks, minor=True)
			ax.tick_params(axis='both', which='minor', labelsize=6)
	fig.canvas.draw()


def get_plot_args(algo_list, hue, plot_args, plot_type, x, y):
	if plot_type == PlotType.line:
		default_plot_args = {
			'markers': True,
			'linewidth': 2,
			'ci': None,
			'style': hue,
		}
		if len(algo_list) > 4:
			default_plot_args['dashes'] = False
		if len(algo_list) > 8:
			default_plot_args['markers'] = False
	elif plot_type == PlotType.swarm:
		default_plot_args = {
			'size': 2,
		}
	elif plot_type is PlotType.bar:
		default_plot_args = {
		}
	elif plot_type is PlotType.violin:
		default_plot_args = {
			'scale': 'count',
			'linewidth': 1,
			'inner': None,
			'cut': 0,
		}
	else:
		raise RuntimeError('No valid plot type!')

	plot_args = plot_args or {}
	plot_args = {
		**default_plot_args,
		'x': x,
		'y': y,
		'hue': hue,
		**plot_args,
	}
	if plot_args['hue'] == 'algo':
		plot_args = {'hue_order': algo_list, **plot_args}
		if plot_type is PlotType.line:
			plot_args = {'style_order': algo_list, **plot_args}
		elif plot_args['x'] == 'algo' and plot_type is PlotType.swarm:
			plot_args = {'order': algo_list, **plot_args}
	return plot_args


def get_algo_list(algo_matches, add_algos, data):
	algo_set = set()
	if algo_matches:
		assert (isinstance(algo_matches, list))
		for algo_match in algo_matches:
			algo_data = data.query('algo.str.contains(@algo_match)')
			algo_set = algo_set.union(set(get_unique_values(algo_data, 'algo')))
	else:
		algo_set = set(get_unique_values(data, 'algo'))
	if add_algos:
		for algo in add_algos:
			algo_set.add(algo)
	# Set ground truth as last algorithm
	gt = 'Ground_Truth'
	if gt in algo_set:
		algo_set.remove(gt)
		algo_list = sorted(list(algo_set)) + [gt]
	else:
		algo_list = sorted(list(algo_set))
	return algo_list


def set_xticklabels(ax, xlabel):
	labels = ax.get_xticklabels()

	def cut_label(label, xlabel):
		begin = label.find(xlabel) + len(xlabel)
		end = label.find('_', begin)
		if end == -1:
			return label[begin:]
		return label[begin:end]

	labels = [cut_label(l.get_text(), xlabel) for l in labels]
	ax.set_xticklabels(labels)
