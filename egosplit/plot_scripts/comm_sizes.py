import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import *


def comm_sizes_lfr(data):
	comm_sizes_filter(data, 'LFR_om')


def comm_sizes_filter(data, graphs):
	graph_names = sorted(list(set(data["cover_comm_sizes"].query("graph.str.contains(@graphs)")["graph"])))
	algos = ["ground_truth"] + algo_sets["ego_parameters"]

	for graph_name in graph_names:
		filtered_data = data["cover_comm_sizes"].query(("graph == @graph_name and algo in @algos"))
		fig, ax = plt.subplots()
		# sns.violinplot(x="graph",
		# 			   y="comm_size",
		# 			   hue="algo",
		# 			   scale="count",
		# 			   linewidth=0.5,
		# 			   inner=None,
		# 			   palette="bright",
		# 			   data=filtered_data,
		# 			   ax=ax,
		# 			   )
		sns.swarmplot(
			x="algo",
			y="comm_size",
			hue="algo",
			# dodge=True,
			order=algos,
			hue_order=algos,
			# palette=colors,
			size=2,
			data=filtered_data,
			ax=ax,
		)
		sns.despine(ax=ax)
		ax.set(ylabel="size (log2)",
		       ylim=(0,10))
		ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
		fig.suptitle("Size of communities, " + graph_name)
		plt.tight_layout(rect=(0, 0, 1, 0.96))
		fig.savefig(file_prefix + "communities/" + 'comm_sizes_' + graph_name + '.pdf')
	# plt.close(fig)


def comm_sizes_ego(data):
	comm_sizes_ego_filter(data, 'LFR_om')
	comm_sizes_ego_filter(data, 'LFR_mu')


def comm_sizes_ego_filter(data, graphs):
	filtered_data = data["cover_comm_sizes"].query(
		"(algo == 'EgoSplitting_(PLM)' or algo == 'ground_truth')"
		" and graph.str.contains(@graphs)")
	fig, ax = plt.subplots()
	sns.violinplot(x="graph",
	               y="comm_size",
	               hue="algo",
	               split=True,
	               scale="count",
	               linewidth=1,
	               inner=None,
	               palette="bright",
	               data=filtered_data,
	               ax=ax)
	sns.despine(ax=ax)
	ax.set(
		ylabel="size (log2)",
		xticklabels=[label.get_text()[4:] for label in ax.get_xticklabels()],
	)
	ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
	fig.suptitle("Size of coummunities, LFR graphs")
	plt.tight_layout(rect=(0, 0, 1, 0.96))
	fig.savefig(file_prefix + "communities/" + 'comm_sizes_' + graphs + '_ego.pdf')
# plt.close(fig)
