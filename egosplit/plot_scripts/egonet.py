import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import *
from .utils import filter_data


def create_egonet_plots(data):
	# egonet_comm_partition(data, "LFR_om_5", "Ego_PLM_1.0_", "")
	# egonet_partition_composition(data, "LFR_om_5", "Ego_PLM_1.0_", "")
	metrics = [
		"community_cohesion",
		"partition_exclusivity",
		# "ego_partitioning_score",
		# "merged_external_nodes",
		# "parts_per_comm",
		# "comms_per_part",
	]
	egonet_metrics(data, graphs="LFR_om", xlabel="om", algos="PLM_1.0_",
	               remove_algo_part="Ego_PLM_1.0_", file_name_add="", title_add=", PLM_1.0",
	               metrics=metrics)
	egonet_metrics(data, graphs="LFR_om", xlabel="om", algos="triangles",
               remove_algo_part="_triangles", file_name_add="_triangles_comm", title_add=", triangles",
               metrics=["community_cohesion"])
	egonet_metrics(data, graphs="LFR_om", xlabel="om", algos="triangles",
	               remove_algo_part="_triangles", file_name_add="_triangles_part", title_add=", triangles",
	               metrics=["partition_exclusivity"])
	egonet_metrics(data, graphs="LFR_om", xlabel="om", algos="triangles",
	               remove_algo_part="_triangles", file_name_add="_triangles_score", title_add=", triangles",
	               metrics=["ego_partitioning_score"])
	egonet_metrics(data, graphs="LFR_om", xlabel="om", algos="triangles",
	               remove_algo_part="_triangles", file_name_add="_triangles_ppc", title_add=", triangles",
	               metrics=["parts_per_comm"])
	egonet_metrics(data, graphs="LFR_om", xlabel="om", algos="triangles",
	               remove_algo_part="_triangles", file_name_add="_triangles_cpp", title_add=", triangles",
	               metrics=["comms_per_part"])
	egonet_metrics(data, graphs="LFR_om", xlabel="om", algos="triangles",
	               remove_algo_part="_triangles", file_name_add="_triangles_ext", title_add=", triangles",
	               metrics=["merged_external_nodes"])
	pass


def egonet_metrics(data, graphs, xlabel, algos, remove_algo_part, file_name_add, title_add, metrics):
	filtered_data = filter_data(data["ego_net_metrics"], graphs, algos)
	if len(filtered_data) is 0:
		return
	filtered_data.query("metric_name in @metrics", inplace=True)

	# if type(algos) == list:
	# 	algo_list = algos
	# else:
	# 	algo_list = sorted(list(set(filtered_data["algo"])))

	fig, ax = plt.subplots()
	sns.lineplot(
		x="graph",
		y="value",
		hue="algo",
		style="metric_name",
		markers=True,
		linewidth=2,
		# palette=colors,
		ci=None,
		data=filtered_data,
		ax=ax,
	)

	sns.despine(ax=ax)
	ax.set(
		xlabel=xlabel,
		ylim=(0, 1.05),
	)
	set_xticklabels(ax, xlabel)
	fig.suptitle("Ego-Metrics, " + graphs + title_add)
	legend_handles, legend_labels = ax.get_legend_handles_labels()
	if type(algos) == str and remove_algo_part == "":
		remove_algo_part = algos
	legend_labels = [l.replace(remove_algo_part, '') for l in legend_labels]
	set_layout(ax, legend_handles, legend_labels)
	file_name = file_prefix + "ego_partition/" + "metrics_" + graphs + file_name_add
	fig.savefig(file_name + ".pdf")

# TODO: ein Wert für Partitions/Community und Communities/Partition
def egonet_comm_partition(data, graphs, algos, name):
	filtered_data = filter_data(data["ego_net_communities"], graphs, algos)
	if len(filtered_data) is 0:
		return

	# if match_algo_name:
	# 	algo_list = sorted(list(set(filtered_data["algo"])))

	# Nodes in other partitions (absolute)
	fig, ax = plt.subplots()
	sns.lineplot(
		x="comm_size",
		y="wrong_nodes",
		hue="algo",
		linewidth=2,
		# palette=colors,
		ci=None,
		data=filtered_data,
		ax=ax,
	)
	sns.despine(left=True, ax=ax)
	plt.suptitle("Number of nodes in other partitions (for each community)\n Graph: " + graphs)
	legend_handles, legend_labels = ax.get_legend_handles_labels()
	if match_algo_name:
		legend_labels = [l.replace(algos, '') for l in legend_labels]
	set_layout(ax, legend_handles, legend_labels)
	fig.savefig(file_prefix + 'ego_partition/' + 'comm_' + graphs + '_total.pdf')

	# Nodes in other partitions (percentage)
	fig, ax = plt.subplots()
	sns.lineplot(
		x="comm_size",
		y="wrong_percentage",
		hue="algo",
		linewidth=2,
		# palette=colors,
		ci=None,
		data=filtered_data,
		ax=ax,
	)
	sns.despine(left=True, ax=ax)
	plt.suptitle("Percentage of nodes in other partitions (for each community)\n Graph: " + graphs)
	legend_handles, legend_labels = ax.get_legend_handles_labels()
	if match_algo_name:
		legend_labels = [l.replace(algos, '') for l in legend_labels]
	set_layout(ax, legend_handles, legend_labels)
	fig.savefig(file_prefix + 'ego_partition/' + 'comm_' + graphs + '_percent.pdf')

	# Number of partitions
	# TODO: Swarm Plot, für alle comm_sizes (Größe egal)
	fig, ax = plt.subplots()
	sns.lineplot(
		x="comm_size",
		y="num_partitions",
		hue="algo",
		linewidth=2,
		# palette=colors,
		ci=None,
		data=filtered_data,
		ax=ax,
	)
	sns.despine(left=True, ax=ax)
	plt.suptitle("Number of partitions in which the community is dominant\n Graph: " + graphs)
	legend_handles, legend_labels = ax.get_legend_handles_labels()
	if match_algo_name:
		legend_labels = [l.replace(algos, '') for l in legend_labels]
	set_layout(ax, legend_handles, legend_labels)
	fig.savefig(file_prefix + 'ego_partition/' + 'comm_' + graphs + '_partcnt.pdf')

	# # Kernel density
	# for algo in filtered_data.groupby("algo").mean().index.values:
	# 	grid = sns.jointplot(
	# 		x="comm_size",
	# 		y="wrong_percentage",
	# 		kind="kde",
	# 		data=filtered_data.query("algo == @algo"),
	# 		xlim=(0, 100),
	# 		ylim=(-0.1, 1),
	# 	)
	# 	plt.suptitle(algo)
	# 	plt.suptitle("Percentage of nodes in other partitions (for each community)"
	# 	             "\n Graph: " + graphs + ", Algo: " + algo)
	# 	# set_layout(ax)
	# 	grid.savefig(file_prefix + 'ego_partition/' + 'comm_' + graphs + '_partcnt' + algo
	# 	            + '.pdf')

	# TODO: Summe comm_size / Summe wrong_nodes -> total wrong percentage
	# Das gleich bei Partition

def egonet_partition_composition(data, graphs, algos, name):
	filtered_data = filter_data(data["ego_net_partitions"], graphs, algos,)
	if len(filtered_data) is 0:
		return

	# Wrong nodes in partition
	for y_val in ["wrong_nodes", "wrong_percentage", "wrong_percentage_gt", "wrong_percentage_other"]:
		fig, ax = plt.subplots()
		sns.lineplot(
			x="partition_size",
			y=y_val,
			hue="algo",
			# style="algo",
			# markers=True,
			linewidth=1.5,
			# palette=colors,
			ci=None,
			data=filtered_data,
			ax=ax,
		)
		sns.despine(left=True, ax=ax)
		plt.suptitle("Wrong nodes in the partition, Graph: " + graphs)
		legend_handles, legend_labels = ax.get_legend_handles_labels()
		if match_algo_name:
			legend_labels = [l.replace(algos, '') for l in legend_labels]
		set_layout(ax, legend_handles, legend_labels)
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
	# 	legend_handles, legend_labels = ax.get_legend_handles_labels()
	# 	if match_algo_name:
	# 		legend_labels = [l.replace(algos, '') for l in legend_labels]
	# 	set_layout(ax, legend_handles, legend_labels)
	# 	plt.savefig(file_prefix + 'ego_partition/' + 'part_' + 'kde_' + algo + '.pdf')

	for ylim in [True, False]:
		fig, ax = plt.subplots()
		for algo_name in filtered_data.groupby("algo").mean().index.values:
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
			ax.set(ylim=(0, 160))
		plt.suptitle("Partition sizes, Graph: " + graphs)
		legend_handles, legend_labels = ax.get_legend_handles_labels()
		if match_algo_name:
			legend_labels = [l.replace(algos, '') for l in legend_labels]
		set_layout(ax, legend_handles, legend_labels)
		plt.savefig(file_prefix + 'ego_partition/' + 'part_sizes_' + str(int(ylim)) + '.pdf')
