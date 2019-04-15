import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import *


def partition_counts(data):
	partition_counts_filter(data, 'LFR_om', 'om', "", "", True)


# partition_counts_filter('LFR_mu')


def partition_counts_filter(data, graphs, xlabel, algos, name, match_algo_name=False):
	filtered_data = data["ego_net_partition_counts"].query("graph.str.contains(@graphs)")
	if match_algo_name:
		filtered_data = filtered_data.query("algo.str.contains(@algos)")
		algos = filtered_data.groupby("algo").mean().index.values
	else:
		filtered_data = filtered_data.query("algo in @algos")

	for algo in algos:
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
		             data=filtered_data.query("algo == @algo"),
		             ax=ax)
		sns.despine(left=True, ax=ax)
		ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3,
		          prop={'size': 9})
		fig.suptitle("Partition counts per node (avg)")
		plt.tight_layout(rect=(0, 0, 1, 0.96))
		fig.savefig(
			file_prefix + 'partition_counts/' + graphs + "_" + algo + '.pdf')
# plt.close(fig)
