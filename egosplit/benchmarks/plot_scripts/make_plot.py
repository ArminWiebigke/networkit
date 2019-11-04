import os
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns

from enum import Enum
from matplotlib.ticker import MultipleLocator, StrMethodFormatter, Formatter
from pandas import DataFrame, Series

from egosplit.benchmarks.plot_scripts.plot_config import set_layout, set_legend, legend_font_size
from egosplit.benchmarks.plot_scripts.read_data import create_column_if_missing


class PlotType(Enum):
	line = 1
	bar = 2
	swarm = 3
	violin = 4


# Create a plot
def make_plot(data,
              output_dir,
              plot_subdir,
              file_name,
              legend_file_name,
              x,
              y,
              hue,
              plot_args,
              ax_set,
              filter_data=None,
              x_filter=None,
              graph_filter='',
              algo_matches='',
              replace_legend=None,
              add_algos=None,
              remove_algo_part=None,
              one_plot_per_graph=False,
              use_graph_id=False,
              title=None,
              max_legend_cols=5,
              max_legend_width=90,
              plot_type=PlotType.line,
              legend_in_separate_file=False):
	replace_legend = replace_legend or {}
	add_algos = add_algos or []
	remove_algo_part = remove_algo_part or []

	# Create the columns for x, y, hue if not already present in data
	for column in [x, y, hue]:
		create_column_if_missing(data, column)
	# Filter data
	if isinstance(graph_filter, str):
		filtered_data = data.query('`Graph Name`.str.contains(@graph_filter)').copy()
	else:
		assert (isinstance(graph_filter, list))
		filtered_data = data.query('`Graph Name` in @graph_filter').copy()
	if filter_data:
		filtered_data.query(filter_data, inplace=True)
	if x_filter:
		filtered_data.query(x_filter, inplace=True)
	graphs = get_unique_values(filtered_data, 'Graph Name')
	filter_algos = [a[1:] for a in algo_matches if len(a) > 0 and a[0] == '!']
	algo_matches = [a for a in algo_matches if len(a) == 0 or a[0] != '!']
	algo_list = get_algo_list(algo_matches, add_algos, filter_algos, filtered_data)

	filtered_data.query('`Algorithm` in @algo_list', inplace=True)
	if len(filtered_data) is 0:
		return

	# Create plots
	def create_plot(graph_data, remove_xlabel=False):
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
		set_ax(fig, ax, ax_set, x)

		def replace_label_func(l):
			return replace_legend_entry(l, remove_algo_part, replace_legend)
		handles, labels = ax.get_legend_handles_labels()
		remove_label(hue, handles, labels)
		labels = [replace_label_func(l) for l in labels]
		handles, labels, num_columns = transpose_legend(handles, labels, max_legend_cols, max_legend_width)
		set_legend(ax, num_columns, handles, labels)

		if legend_in_separate_file:
			ax.legend().remove()
			legend_file = output_dir + plot_subdir + '{}.pdf'.format(legend_file_name)
			save_legend(handles, labels, num_columns, legend_file)
		if remove_xlabel:
			ax.get_xaxis().set_visible(False)
		set_layout()

		if title:
			fig.suptitle(title)
		plot_file = output_dir + plot_subdir + file_name + '.pdf'
		fig.savefig(plot_file)
		# crop_pdf(plot_file)
		plt.close(fig)

		with open(output_dir + 'tables/' + plot_subdir + file_name + '.tex', 'w') as table_file:
			table_file.write(latex_table(graph_data, hue, x, y, replace_label_func))

	# Make one plot per graph if graph is not on the x-axis
	if one_plot_per_graph:
		file_name_base = file_name
		for graph in graphs:
			file_name = file_name_base + '_' + graph
			graph_data = filtered_data.query('`Graph Name` == @graph')
			if use_graph_id:
				graph_ids = get_unique_values(graph_data.query('Algorithm != \'Ground Truth\''),
				                              'Graph ID')
				if len(graph_ids) == 0:
					continue
				first_graph_id = graph_ids[0]
				graph_data = graph_data.query('`Graph ID` == @first_graph_id')
			remove_xlabel = x == 'Graph Name'
			create_plot(graph_data, remove_xlabel)
	else:
		create_plot(filtered_data)


def remove_label(to_remove, handles, labels):
	if to_remove in labels:
		idx = labels.index(to_remove)
		del labels[idx]
		del handles[idx]


def save_legend(handles, labels, num_columns, legend_file):
	if os.path.exists(legend_file):
		return
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.legend(
		handles=handles,
		labels=labels,
		ncol=num_columns,
		fontsize=legend_font_size(),
	)
	# hide the axes frame and the x/y labels
	ax.axis('off')
	fig.savefig(legend_file)
	crop_pdf(legend_file)
	plt.close(fig)


def crop_pdf(file_name):
	subprocess.call(["pdfcrop",
	                 "--margins", "0 0 0 0",
	                 file_name, file_name])


def transpose_legend(handles, labels, max_columns, max_legend_width):
	num_columns = max_columns
	while legend_too_long(labels, num_columns, max_legend_width):
		num_columns -= 1
	new_handles = []
	new_labels = []
	for i in range(num_columns):
		new_handles.extend(handles[i::num_columns])
		new_labels.extend(labels[i::num_columns])
	return new_handles, new_labels, num_columns


def legend_too_long(labels, num_columns, max_width):
	if num_columns == 1:
		return False
	num_columns = min(num_columns, len(labels))
	max_width = max_width - 9 * num_columns  # Extra space between columns
	column_length_sum = 0
	for i in range(num_columns):
		labels_in_column = [labels[idx] for idx in range(i, len(labels), num_columns)]
		column_length_sum += max([len(l) for l in labels_in_column])
	if column_length_sum > max_width:
		return True
	return False


def replace_legend_entry(label, remove_algo_part, replace_legend):
	for remove in remove_algo_part:
		label = label.replace(remove, '')
	for remove, with_value in replace_legend.items():
		if label.strip(' ') == remove:
			label = with_value
			break
	# label = label.replace('_', ' ')
	# label = label.replace('&', '+')
	return label


def latex_table(data, hue, x, y, replace_label_func):
	hue_values = sorted(get_unique_values(data, hue))
	x_values = sorted(get_unique_values(data, x))
	table_data = DataFrame(columns=x_values)
	for hue_value in hue_values:
		hue_data = data.query('`{}` == @hue_value'.format(hue))
		mean_data = hue_data.groupby(x).mean()
		name = hue_value
		if hue == 'Algorithm':
			name = replace_label_func(name)
		series = Series(mean_data[y], name=name)
		table_data = table_data.append(series)
	x_values = [str(remove_facebook_prefix(l)) for l in x_values]
	return table_data.to_latex(float_format=float_format, escape=False, header=x_values)


def float_format(x, decimals=3):
	format_str = '{:.' + str(decimals) + 'f}'
	return format_str.format(x)


def confirm_plot_type(plot_type, graph_data, x_values):
	if plot_type == PlotType.swarm and len(graph_data) > 2000:
		return PlotType.violin

	if x_values == 1 and plot_type == PlotType.line:
		plot_type = PlotType.bar
	return plot_type


def get_unique_values(filtered_data, column):
	return filtered_data.groupby(column).mean().index.values


def draw_plot(plot_type, plot_args):
	plot_functions = {
		PlotType.line: sns.lineplot,
		PlotType.bar: sns.barplot,
		PlotType.swarm: sns.swarmplot,
		PlotType.violin: sns.violinplot,
	}
	return plot_functions[plot_type](**plot_args)


class MinorFormatter(Formatter):
	def __call__(self, x, pos=None):
		s = '{:.1f}'.format(x)
		return s[-2:]


def set_ax(fig, ax, ax_set, x):
	ax_set = ax_set or {}
	ax.set(
		**ax_set,
	)
	if x == 'Communities per Node' and ax.lines:
		ax.xaxis.set(
			major_locator=MultipleLocator(1),
			major_formatter=StrMethodFormatter('{x:.0f}'),
			minor_formatter=MinorFormatter(),
		)
		minor_ticks = [x_val for x_val in ax.lines[0].get_xdata() if 1 < x_val < 2]
		if minor_ticks:
			ax.set_xticks(minor_ticks, minor=True)
			ax.tick_params(axis='both', which='minor', labelsize=9)

	# for barplot
	xlabels = [l.get_text() for l in ax.get_xticklabels()]
	if 'FB_' in xlabels[0]:
		xlabels = [remove_facebook_prefix(l) for l in xlabels]
	if xlabels[0] != '':
		xlabels = [l.replace("_", ' ') for l in xlabels]
		ax.set_xticklabels(xlabels)
	# for lineplot
	if ax.lines:
		ticks = [x_val for x_val in ax.lines[0].get_xdata()]
		if isinstance(ticks[0], str) and 'FB_' in ticks[0]:
			ax.set_xticks(ax.get_xticks())
			ax.set_xticklabels([remove_facebook_prefix(t) for t in ticks])
	fig.canvas.draw()


def remove_facebook_prefix(label):
	if not isinstance(label, str):
		return label
	if 'FB_' in label:
		label = label[5:]
	return label


def get_plot_args(algo_list, hue, plot_args, plot_type, x, y):
	if plot_type == PlotType.line:
		default_plot_args = {
			'markers': True,
			'markersize': 9,
			'linewidth': 2,
			'ci': None,
			'style': hue,
			'dashes': False,
		}
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
	if plot_args['hue'] == 'Algorithm':
		plot_args = {'hue_order': algo_list, **plot_args}
		if plot_type is PlotType.line:
			plot_args = {'style_order': algo_list, **plot_args}
		elif plot_args['x'] == 'Algorithm' and plot_type is PlotType.swarm:
			plot_args = {'order': algo_list, **plot_args}
	return plot_args


def get_algo_list(algo_matches, add_algos, filter_algos, data):
	algo_set = set()
	assert (isinstance(algo_matches, list))
	for algo_match in algo_matches:
		algo_match = algo_match.replace('+', r'\+')  # Pandas uses regex syntax in str.contains
		algo_data = data.query('Algorithm.str.contains(@algo_match)')
		algo_set = algo_set.union(set(get_unique_values(algo_data, 'Algorithm')))
	# algo_set = set(get_unique_values(data, 'Algorithm'))
	for filter in filter_algos:
		algo_set = set(a for a in algo_set if filter not in a)
	for algo in add_algos:
		algo_set.add(algo)
	# Set ground truth as last algorithm
	gt = 'Ground Truth'
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
