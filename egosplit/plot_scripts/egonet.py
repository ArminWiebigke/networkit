import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import *


def egonet_comm_partition(data, graphs):
	algos = ["ground_truth"] + algo_sets["ego_parameters"]
	filtered_data = data["ego_net_partitions"].query("graph.str.contains(@graphs)"
	                                                 " and algo in @algos")

	# Nodes in other partitions (absolute)
	fig, ax = plt.subplots()
	sns.lineplot(
		x="comm_size",
		y="wrong_nodes",
		hue="algo",
		linewidth=3,
		# palette=colors,
		ci=None,
		data=filtered_data,
		ax=ax,
	)
	sns.despine(left=True, ax=ax)
	ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
	plt.suptitle("Number of nodes in other partitions (for each community)\n Graph: " + graphs)
	plt.tight_layout(rect=(0, 0, 1, 0.96))
	fig.savefig(file_prefix + 'ego_partition/' + 'comm_' + graphs + '_total.pdf')

	# Nodes in other partitions (percentage)
	fig, ax = plt.subplots()
	sns.lineplot(
		x="comm_size",
		y="wrong_percentage",
		hue="algo",
		linewidth=3,
		# palette=colors,
		ci=None,
		data=filtered_data,
		ax=ax,
	)
	sns.despine(left=True, ax=ax)
	ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
	plt.suptitle("Percentage of nodes in other partitions (for each community)\n Graph: " + graphs)
	plt.tight_layout(rect=(0, 0, 1, 0.96))
	fig.savefig(file_prefix + 'ego_partition/' + 'comm_' + graphs + '_percent.pdf')

	# Number of partitions
	fig, ax = plt.subplots()
	sns.lineplot(
		x="comm_size",
		y="num_partitions",
		hue="algo",
		linewidth=3,
		# palette=colors,
		ci=None,
		data=filtered_data,
		ax=ax,
	)
	sns.despine(left=True, ax=ax)
	ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
	plt.suptitle("Number of partitions in which the community is dominant\n Graph: " + graphs)
	plt.tight_layout(rect=(0, 0, 1, 0.96))
	fig.savefig(file_prefix + 'ego_partition/' + 'comm_' + graphs + '_partcnt.pdf')

# # Kernel density
# for algo in filtered_data.groupby("algo").mean().index.values:
# 	sns.jointplot(
# 		x="comm_size",
# 		y="wrong_nodes",
# 		kind="kde",
# 		data=filtered_data.query("algo == @algo"),
# 		xlim=(0, 38),
# 		ylim=(-1, 20),
# 	)
# 	plt.suptitle(algo)
# 	plt.tight_layout(rect=(0, 0, 1, 0.96))


def egonet_partition_composition(data, graphs):
	algos = ["ground_truth"] + algo_sets["ego_parameters"]
	filtered_data = data["ego_net_partition_composition"].query("graph.str.contains(@graphs)"
	                                                            " and algo in @algos")

	# Wrong nodes in partition
	for y_val in ["wrong_nodes", "wrong_percentage", "wrong_percentage_gt", "wrong_percentage_other"]:
		fig, ax = plt.subplots()
		sns.lineplot(
			x="partition_size",
			y=y_val,
			hue="algo",
			# style="algo",
			# markers=True,
			linewidth=2.5,
			# palette=colors,
			ci=None,
			data=filtered_data,
			ax=ax,
		)
		sns.despine(left=True, ax=ax)
		ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
		plt.suptitle("Wrong nodes in the partition, Graph: " + graphs)
		plt.tight_layout(rect=(0, 0, 1, 0.96))
		fig.savefig(file_prefix + 'ego_partition/' + 'part_' + graphs + '_' + y_val + '.pdf')

	# # Partition sizes
	# fig, ax = plt.subplots()
	# sns.violinplot(x="graph",
	#                y="partition_size",
	#                hue="algo",
	#                scale="count",
	#                linewidth=1,
	#                inner=None,
	#                palette=colors,
	#                data=filtered_data,
	#                ax=ax)

	# # Kernel density
	# for algo in filtered_data.groupby("algo").mean().index.values:
	# 	sns.jointplot(
	# 		x="partition_size",
	# 		y="wrong_nodes",
	# 		kind="kde",
	# 		data=filtered_data.query("algo == @algo"),
	# 		xlim=(-16, 80),
	# 		ylim=(-10, 50),
	# 	)
	# 	plt.suptitle("Number of wrong nodes in the partition, Graph: " + graphs + ", "+ algo)
	# 	plt.tight_layout(rect=(0, 0, 1, 0.96))
	# 	plt.savefig(file_prefix + 'ego_partition/' + 'part_' + 'kde_' + algo + '.pdf')

	for ylim in [True, False]:
		fig, ax = plt.subplots()
		for algo_name in data["ego_net_partition_composition"].groupby("algo").mean().index.values:
			count_communities = filtered_data.query("algo == @algo_name")
			sns.distplot(
				a=count_communities["partition_size"],
				# hist=False,
				kde=False,
				bins=range(0, 60, 1),
				hist_kws={"histtype": "step",
				          "linewidth": 2,
				          "alpha": 1,
				          #"color": colors[algo_name],
				          },
				label=algo_name,
				ax=ax,
			)
		sns.despine(left=True, ax=ax)
		if ylim:
			ax.set(ylim=(0, 900))
		ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
		plt.suptitle("Partition sizes, Graph: " + graphs)
		plt.tight_layout(rect=(0, 0, 1, 0.96))
		plt.savefig(file_prefix + 'ego_partition/' + 'part_sizes_' + str(int(ylim)) + '.pdf')
