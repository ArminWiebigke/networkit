import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import *
from .utils import filter_data


def create_comm_sizes_plots(data):
	algos = [
		"ground_truth",
		"Ego_PLM_1.0_base",
		"Ego_PLM_1.0_edges",
		"Ego_PLM_1.0_triangles_minTri0",
		"Ego_PLM_1.0_triangles_minTri2",
		"Ego_PLM_1.0_base_clean_OSLOM",
		"Ego_PLM_1.0_edges_clean_OSLOM",
		"Ego_PLM_1.0_triangles_minTri0_clean_OSLOM",
		"Ego_PLM_1.0_triangles_minTri2_clean_OSLOM",
	]
	comm_sizes_filter(data, 'LFR_om', '', algos, "Ego_PLM_1.0_", "",)
	pass


def comm_sizes_filter(data, graphs, xlabel, algos, remove_algo_part, name):
	filtered_data = filter_data(data["cover_comm_sizes"], graphs, algos)
	if len(filtered_data) is 0:
		return
	graph_names = sorted(list(set(filtered_data["graph"])))
	if type(algos) == list:
		algo_list = algos
	else:
		algo_list = sorted(list(set(filtered_data["algo"])))

	for graph_name in graph_names:
		graph_data = filtered_data.query("graph == @graph_name")
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
			order=algo_list,
			hue_order=algo_list,
			# palette=colors,
			size=2,
			data=graph_data,
			ax=ax,
		)
		sns.despine(ax=ax)
		ax.set(
			xlabel=xlabel,
			ylabel="size (log2)",
			ylim=(0, 10),
			xticklabels=[""]
		)

		fig.suptitle("Size of communities, " + graph_name)
		legend_handles, legend_labels = ax.get_legend_handles_labels()
		if type(algos) == str and remove_algo_part == "":
			remove_algo_part = algos
		legend_labels = [l.replace(remove_algo_part, '') for l in legend_labels]
		set_layout(ax, legend_handles, legend_labels)

		file_name = file_prefix + "communities/" + 'comm_sizes_' + graph_name + name
		fig.savefig(file_name + ".pdf")
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
