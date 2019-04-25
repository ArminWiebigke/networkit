import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import *
from .utils import filter_data


def create_partition_counts_plots(data):
	# partition_counts_per_algo(data, 'LFR_om', 'om', "", "", True)
	# partition_counts_all_algos(data, 'LFR_om', 'om', "_triangles", "")
	pass



def partition_counts_per_algo(data, graphs, xlabel, algos, name):
	filtered_data = filter_data(data["ego_net_metrics"], graphs, algos)
	if len(filtered_data) is 0:
		return
	if type(algos) == list:
		algo_list = algos
	else:
		algo_list = sorted(list(set(filtered_data["algo"])))

	for algo in algo_list:
		algo_data = filtered_data.query("algo == @algo")
		fig, ax = plt.subplots()
		sns.lineplot(x="graph",
		             y="count",
		             hue="count_name",
		             # style="algo",
		             # markers=markers,
		             markers=True,
		             dashes=False,
		             linewidth=2,
		             markersize=8,
		             # palette=colors,
		             data=algo_data,
		             ax=ax)
		sns.despine(left=True, ax=ax)
		ax.set(
			xlabel=xlabel,
		)

		set_xticklabels(ax, xlabel)
		fig.suptitle("Partition counts per node (avg), " + algo)
		set_layout(ax)

		file_name = file_prefix + 'partition_counts/' + graphs + "_" + algo + name
		fig.savefig(file_name + ".pdf")
		# plt.close(fig)


def partition_counts_all_algos(data, graphs, xlabel, algos, name):
	filtered_data = filter_data(data["ego_net_ego_metrics"], graphs, algos)
	metric_name = "num_three_plus_partitions"
	filtered_data = filtered_data.query("metric_name == @metric_name")
	if len(filtered_data) is 0:
		return
	# if match_algo_name:
	# 	algos = filtered_data.groupby("algo").mean().index.values

	fig, ax = plt.subplots()
	sns.lineplot(x="graph",
	             y="value",
	             hue="algo",
	             style="algo",
	             # markers=markers,
	             ci=None,
	             markers=True,
	             dashes=False,
	             linewidth=2,
	             # markersize=8,
	             # palette=colors,
	             data=filtered_data,
	             ax=ax)
	sns.despine(left=True, ax=ax)
	ax.set(
		xlabel=xlabel,
	)

	set_xticklabels(ax, xlabel)
	fig.suptitle("Partition counts per node (avg), " + metric_name)
	set_layout(ax)

	file_name = file_prefix + 'partition_counts/' + graphs + "_" + metric_name + name
	fig.savefig(file_name + ".pdf")
