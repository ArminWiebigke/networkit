import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import *
from .utils import filter_data


def create_metric_plots(data):
	# metrics_filter(data, 'LFR_om', 'om', "Ego_PLM_1.0", "", True)
	metrics_filter(data, 'LFR_om', 'om', "_triangles", "", True)
	# metrics_2dim(data, 'LFR_om', '', '', '', True)
	pass

# TODO: Remove match_algo_name parameter
def metrics_filter(data, graphs, xlabel, algos, name, match_algo_name=False):
	filtered_data = filter_data(data["metrics"], graphs, algos)
	if len(filtered_data) is 0:
		return
	metrics = [
		'time',
		# 'f1',
		# 'f1_rev',
		'nmi',
		# 'entropy',
	]
	for metric in metrics:
		fig, ax = plt.subplots()
		sns.lineplot(x="graph",
		             y=metric,
		             hue="algo",
		             ci=None,
		             style="algo",
		             markers=True,
		             # markers=markers,
		             dashes=False,
		             linewidth=2,
		             markersize=6,
		             # palette=sns.light_palette("green", 20),
		             # palette=colors,
		             data=filtered_data,
		             ax=ax)
		sns.despine(ax=ax)
		ax.set(
			xlabel=xlabel,
			ylabel=metric_names[metric]["y_val"],
		)

		set_xticklabels(ax, xlabel)
		if metric != "time":
			ax.set(ylim=(-0.01, 1))
		fig.suptitle(metric_names[metric]["description"] + ", LFR graphs")
		legend_handles, legend_labels = ax.get_legend_handles_labels()
		if match_algo_name:
			legend_labels = [l.replace(algos, '') for l in legend_labels]
		set_layout(ax, legend_handles, legend_labels)
		file_name = file_prefix + "metrics/" + metric_names[metric]["file_name"] + "_" \
		            + graphs + name
		fig.savefig(file_name + ".pdf")


def metrics_2dim(data, graphs, xlabel, algos, name, match_algo_name=False):
	filtered_data = filter_data(data["metrics"], graphs, algos, match_algo_name)
	if len(filtered_data) is 0:
		return
	filtered_data.query("algo != 'ground_truth'", inplace=True)
	filtered_data.loc[:, 'PLM_gamma'] = filtered_data.loc[:, 'algo'].str[8:11].astype(np.float)
	filtered_data.loc[:, 'triangle_thresh'] = filtered_data.loc[:, 'algo'].str[22:25].astype(np.float)
	if match_algo_name:
		graphs = filtered_data.groupby("graph").mean().index.values

	metrics = [
		# 'time',
		# 'f1',
		# 'f1_rev',
		'nmi',
		# 'entropy',
	]

	#
	hue = "PLM_gamma"
	for graph in graphs:
		for metric in metrics:
			graph_data = filtered_data.query("graph == @graph")
			num_hues = len(filtered_data.groupby(hue).mean().index.values)
			fig, ax = plt.subplots()
			sns.lineplot(x="triangle_thresh",
			             y=metric,
			             hue=hue,
			             ci=None,
			             style=hue,
			             markers=True,
			             linewidth=2,
			             markersize=6,
			             palette=sns.light_palette("green", num_hues),
			             data=graph_data,
			             ax=ax)
			sns.despine(ax=ax)
			ax.set(
				ylabel=metric_names[metric]["y_val"],
				ylim=0,
			)

			fig.suptitle(metric_names[metric]["description"] + ", " + graph)
			legend_handles, legend_labels = ax.get_legend_handles_labels()
			legend_labels = [hue] + [l[:3] for l in legend_labels[1:]]
			set_layout(ax, legend_handles, legend_labels)
			file_name = file_prefix + "metrics/" + "2dim_" + metric_names[metric]["file_name"] + "_" \
			            + graph + name
			fig.savefig(file_name + ".pdf")
