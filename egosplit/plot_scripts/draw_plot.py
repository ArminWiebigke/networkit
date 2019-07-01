import matplotlib.pyplot as plt
import seaborn as sns

from enum import Enum

from .config import file_prefix


class PlotType(Enum):
	line = 1
	bar = 2
	swarm = 3
	violin = 4


# Create a plot
def make_plot(data,
              graph_filter='',
              algo_matches='',
              add_algos=None,
              xlabel=None,
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
	# Filter data
	if isinstance(graph_filter, str):
		filtered_data = data.query('graph.str.contains(@graph_filter)').copy()
		graphs = get_unique_values(filtered_data, 'graph')
	else:
		assert (isinstance(graph_filter, list))
		filtered_data = data.query('graph in @graph_filter').copy()
		graphs = graph_filter
	algo_list = get_algo_list(algo_matches, add_algos, filtered_data)
	filtered_data.query('algo in @algo_list', inplace=True)
	if len(filtered_data) is 0:
		return

	create_new_columns(filtered_data, hue, x)
	if type(x) != str and x:
		x = x['name']
	if type(hue) != str and hue:
		hue = hue['name']

	plot_args = get_plot_args(algo_list, hue, plot_args, plot_type, x, y)

	# Create plots
	def create_plot(graph_data, graph_name):
		fig, ax = plt.subplots()
		this_plot_args = {
			**plot_args,
			'data': graph_data,
			'ax': ax,
		}
		draw_plot(plot_type, this_plot_args)

		sns.despine(ax=ax)
		clean_xlabel(ax, ax_set, fig, xlabel)
		clean_legend(algo_matches, ax, remove_algo_part)
		add_title = ', ' + graph_name if graph_name else ""
		fig.suptitle(title + add_title)
		add_file = '_' + graph_name if graph_name else ""
		fig.savefig(file_prefix + file_name + add_file + '.pdf')
		plt.close(fig)

	# Make one plot per graph if graph is not on the x-axis
	if plot_args['x'] != 'graph' or one_plot_per_graph:
		for graph in graphs:
			graph_data = filtered_data.query('graph == @graph')
			create_plot(graph_data, graph)
	else:
		if (isinstance(graph_filter, str)):
			graph_name = graph_filter
		else:
			graph_name = ""
		create_plot(filtered_data, graph_name)


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


def clean_xlabel(ax, ax_set, fig, xlabel):
	ax_set = ax_set or {}
	old_xlabel = ax.get_xlabel()
	xlabel = xlabel or old_xlabel
	ax_set = {
		'xlabel': xlabel,
		**ax_set,
	}
	ax.set(
		**ax_set,
	)
	fig.canvas.draw()
	if xlabel != old_xlabel:
		set_xticklabels(ax, xlabel + '_')


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


def create_new_columns(graph_filtered_data, hue, x):
	# Create new data columns for x or hue if necessary
	# Example: Create the x column by taking the value between '_f-' and the next '*'
	#   from the 'algo' column
	# x = {'name': 'factor', 'create_from': 'algo', 'str_start': '_f-', 'str_end': '*'}
	new_columns = []
	if type(x) != str and x:
		new_columns.append(x)
	if type(hue) != str and hue:
		new_columns.append(hue)
	for new_column in new_columns:
		for row in graph_filtered_data.itertuples(index=True, name='Pandas'):
			index = row.Index
			from_string = graph_filtered_data.loc[index, new_column['create_from']]
			start_idx = from_string.find(new_column['str_start']) + len(
				new_column['str_start'])
			if new_column['str_end'] == '':
				end_idx = len(from_string)
			else:
				end_idx = from_string.find(new_column['str_end'], start_idx)
			value = from_string[start_idx:end_idx]
			try:
				value = float(value)
			except ValueError:
				pass
			graph_filtered_data.loc[index, new_column['name']] = value


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


def set_layout(ax, legend_handles=None, legend_labels=None):
	legend_args = {
		"loc": "lower center",
		"bbox_to_anchor": (0.5, 1.01),
		"ncol": 3,
		"prop": {'size': 9}
	}
	if legend_handles is None:
		ax.legend(**legend_args)
	else:
		ax.legend(**legend_args, handles=legend_handles, labels=legend_labels)
	plt.tight_layout(rect=(0, 0, 1, 0.96))
