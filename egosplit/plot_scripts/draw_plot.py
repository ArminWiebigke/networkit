import os
import subprocess

import matplotlib.pyplot as plt
import seaborn as sns

from enum import Enum

from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FormatStrFormatter, \
	StrMethodFormatter, Formatter, FuncFormatter
from pandas import DataFrame, Series

from networkit.stopwatch import clockit
from .config import set_layout, set_legend, get_legend_args, set_sns_style
from .read_data import create_column_if_missing


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
              algos_filter='',
              replace_legend=None,
              add_algos=None,
              remove_algo_part=None,
              one_plot_per_graph=False,
              title=None,
              plot_type=PlotType.line,):
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
	algo_list = get_algo_list(algo_matches, add_algos, filtered_data)
	algo_list = [a for a in algo_list if a not in algos_filter]

	filtered_data.query('`Algorithm` in @algo_list', inplace=True)
	if len(filtered_data) is 0:
		return

	# Create plots
	# set_sns_style(font_size)  # TODO: This is not working, does not change font size
	@clockit
	def create_plot(graph_data, single_x_value=False):
		num_x_values = len(get_unique_values(graph_data, x))
		this_plot_type = confirm_plot_type(plot_type, graph_data, num_x_values)
		this_plot_args = get_plot_args(algo_list, hue, plot_args, this_plot_type, x, y)
		a4_dims = (7, 5)
		# fig, ax = plt.subplots(figsize=a4_dims)
		fig, ax = plt.subplots()
		this_plot_args = {
			**this_plot_args,
			'data': graph_data,
			'ax': ax,
		}
		draw_plot(this_plot_type, this_plot_args)

		sns.despine(ax=ax)
		set_ax(fig, ax, ax_set, x)
		# clean_legend(ax, remove_algo_part)
		ax.legend().remove()
		if single_x_value:
			ax.get_xaxis().set_visible(False)
		set_layout()

		if title:
			fig.suptitle(title)
		fig.savefig(output_dir + plot_subdir + file_name + '.pdf')
		plt.close(fig)

		def replace_label_func(l):
			return replace_legend_entry(l, remove_algo_part, replace_legend)
		with open(output_dir + 'tables/' + plot_subdir + file_name + '.tex', 'w') as f:
			f.write(latex_table(graph_data, hue, x, y, replace_label_func))

		legend_file = output_dir + plot_subdir + '{}.pdf'.format(legend_file_name)
		save_legend(ax, legend_file, hue, replace_label_func)

	# Make one plot per graph if graph is not on the x-axis
	if one_plot_per_graph:
		file_name_base = file_name
		# title_base = title or ''
		for graph in graphs:
			file_name = file_name_base + '_' + graph
			graph_data = filtered_data.query('`Graph Name` == @graph')
			remove_xlabel = x == 'Graph Name'
			create_plot(graph_data, remove_xlabel)
	else:
		create_plot(filtered_data)


def save_legend(ax, legend_file, hue, replace_label_func):
	fig_leg = plt.figure()
	ax_leg = fig_leg.add_subplot(111)
	handles, labels = ax.get_legend_handles_labels()
	if hue in labels:
		idx = labels.index(hue)
		del labels[idx]
		del handles[idx]
	# labels, handles = zip(*(sorted(zip(labels, handles), key=lambda t: t)))
	# Remove/replace labels
	labels = [replace_label_func(l) for l in labels]
	num_columns = 5
	while legend_too_long(labels, num_columns):
		num_columns -=1
	handles, labels = transpose_legend(handles, labels, num_columns)
	ax_leg.legend(
		# title=hue,
		handles=handles,
		labels=labels,
		ncol=num_columns,
		**get_legend_args(),
		loc='center')
	# hide the axes frame and the x/y labels
	ax_leg.axis('off')
	fig_leg.savefig(legend_file)
	plt.close(fig_leg)


def transpose_legend(handles, labels, num_columns):
	new_handles = []
	new_labels = []
	for i in range(num_columns):
		new_handles.extend(handles[i::num_columns])
		new_labels.extend(labels[i::num_columns])
	return new_handles, new_labels


def legend_too_long(labels, num_columns):
	if num_columns == 1:
		return False
	max_width = 90 - 9 * num_columns
	for line in [labels[i:i + num_columns] for i in range(0, len(labels), num_columns)]:
		length = sum([len(l) for l in line])
		if length > max_width:
			return True
	return False


def replace_legend_entry(label, remove_algo_part, replace_legend):
	for remove in remove_algo_part:
		label = label.replace(remove, '')
	for remove, with_value in replace_legend.items():
		if label == remove:
			label = with_value
	return label


@clockit
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
	return table_data.to_latex(float_format=float_format, escape=False)


def float_format(x, decimals=3):
	format_str = '{:.' + str(decimals) + 'f}'
	return format_str.format(x)


def confirm_plot_type(plot_type, graph_data, x_values):
	if plot_type == PlotType.swarm and len(graph_data) > 10000:
		return PlotType.violin

	if x_values == 1 and plot_type == PlotType.line:
		plot_type = PlotType.bar
	return plot_type


def get_unique_values(filtered_data, column):
	return filtered_data.groupby(column).mean().index.values


def clean_legend(ax, remove_algo_part):
	legend_handles, legend_labels = ax.get_legend_handles_labels()
	assert (isinstance(remove_algo_part, list))
	remove_list = remove_algo_part
	for remove in remove_list:
		legend_labels = [l.replace(remove, '') for l in legend_labels]
	set_legend(ax, legend_handles, legend_labels)


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
	if x == 'Communities per Node':
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
		ax.set_xticklabels(xlabels)
	# for lineplot
	if ax.lines:
		ticks = [x_val for x_val in ax.lines[0].get_xdata()]
		if isinstance(ticks[0], str) and 'FB_' in ticks[0]:
			ax.set_xticks(ax.get_xticks())
			ax.set_xticklabels([remove_facebook_prefix(t) for t in ticks])

	fig.canvas.draw()


def remove_facebook_prefix(label):
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


def get_algo_list(algo_matches, add_algos, data):
	algo_set = set()
	assert (isinstance(algo_matches, list))
	for algo_match in algo_matches:
		algo_data = data.query('Algorithm.str.contains(@algo_match)')
		algo_set = algo_set.union(set(get_unique_values(algo_data, 'Algorithm')))
	else:
		algo_set = set(get_unique_values(data, 'Algorithm'))
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
