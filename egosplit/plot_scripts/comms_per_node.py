import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import *


def node_comms_lfr(data):
	node_comms(data, "LFR_om")
	node_comms(data, "LFR_mu")


def node_comms(data, graphs):
	graph_names = sorted(list(set(data["cover_node_comms"].query("graph.str.contains(@graphs)")["graph"])))
	for graph_name in graph_names:
		filtered_data = data["cover_node_comms"].query(("graph == @graph_name"))
		fig, ax = plt.subplots()
		sns.violinplot(x="graph",
		               y="num_comms",
		               hue="algo",
		               scale="count",
		               linewidth=1,
		               inner=None,
		               palette="bright",
		               data=filtered_data,
		               ax=ax)
		sns.despine(ax=ax)
		ax.set(ylabel="#communities (log2)")
		ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
		fig.suptitle("Number of communities per node")
		plt.tight_layout(rect=(0, 0, 1, 0.96))
		fig.savefig(file_prefix + "communities/" + 'node_comms_' + graph_name + '.pdf')
	# plt.close(fig)


def node_comms_ego(data):
	node_comms_ego_filter(data, 'LFR_om')
	node_comms_ego_filter(data, 'LFR_mu')


def node_comms_ego_filter(data, graphs):
	filtered_data = data["cover_node_comms"].query(
		"(algo == 'EgoSplitting_(PLM)' or algo == 'ground_truth') and graph.str.contains(@graphs)")
	fig, ax = plt.subplots()
	sns.violinplot(x="graph",
	               y="num_comms",
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
		ylabel="#communities (log2)",
		xticklabels=[label.get_text()[4:] for label in ax.get_xticklabels()],
	)
	ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3)
	fig.suptitle("Number of communities per node, LFR graphs")
	plt.tight_layout(rect=(0, 0, 1, 0.96))
	fig.savefig(file_prefix + "communities/" + 'node_comms_' + graphs + '.pdf')
# plt.close(fig)

