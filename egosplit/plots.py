import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# Read input *******************************************************************
result_dir = "results/"
data = dict()

for name in ["metrics", "execution_info", "ego_net_partitions",
             "ego_net_partition_composition",
             "cover_comm_sizes", "cover_node_comms", "cover_num_comms"]:
	filename = name + '.result'
	try:
		data[name] = pd.read_csv(result_dir + filename, sep="\s+")
	except FileNotFoundError:
		print("File " + filename + " not found")

metrics = ["f1", "f1_rev", "nmi", "time"]
all_graph_names = data["metrics"].groupby("graph").mean().index.values
all_algo_names = data["metrics"].groupby("algo").mean().index.values

# Create Plots *****************************************************************
file_prefix = "plots/"
opacity = 0.7
error_config = {'ecolor': '0.3',
				'capsize': 4}
colors = {
	"ground_truth": "xkcd:black",
	"Ego_PLP": "xkcd:brown",
	"Ego_PLM": "xkcd:green",
	"Ego_PLM_clean_OSLOM": "xkcd:light green",
	"Ego_LPPotts": "xkcd:light blue",
	"Ego_Infomap": "xkcd:red",
	"Ego_Infomap_clean_OSLOM": "xkcd:light red",
	"Ego_Surprise": "xkcd:goldenrod",
	"OLP": "xkcd:magenta",
	"GCE": "xkcd:fuchsia",
	"MOSES": "xkcd:dark teal",
	"OSLOM": "xkcd:blue",
}

markers = {
	"ground_truth": ".",
	"Ego_PLP": "v",
	"Ego_PLM": "s",
	"Ego_PLM_clean_OSLOM": "v",
	"Ego_LPPotts": "*",
	"Ego_Infomap": "X",
	"Ego_Infomap_clean_OSLOM": "*",
	"Ego_Surprise": "o",
	"OLP": "<",
	"GCE": "P",
	"MOSES": "D",
	"OSLOM": "^",
}

metric_names = {
	"f1": {
		"description": "F1 Score",
		"y_val": "F1",
		"file_name": "F1_score"
	},
	"f1_rev": {
		"description": "F1 Score (reversed)",
		"y_val": "F1",
		"file_name": "F1_score_rev"
	},
	"nmi": {
		"description": "NMI Score",
		"y_val": "NMI",
		"file_name": "NMI_score"
	},
	"time": {
		"description": "Running Time",
		"y_val": "time (s)",
		"file_name": "time"
	},
}
sns.set(context="notebook", style="whitegrid", palette="bright", font_scale=0.8)

algo_sets = dict()
algo_sets["clean"] = [
	"Ego_PLM", "Ego_PLM_clean_OSLOM", "OSLOM",
]
algo_sets["ego"] = [
	"Ego_PLP", "Ego_LPPotts", "Ego_PLM", "Ego_Infomap", "Ego_Surprise",
]
algo_sets["comm_sizes"] = [
	"ground_truth", "Ego_PLM", "Ego_Infomap", "Ego_Surprise", "Ego_PLP", "MOSES", "OSLOM", "GCE"
]
algo_sets["ego"] = [
	"Ego_PLM", "Ego_Infomap", "Ego_PLP", "Ego_Surprise"
]


def filter_data(data, conditions):
	for cond in conditions:
		data = data[data[cond[0]] == cond[1]]
	return data


def get_value(base_data, x_name, x_val, y_name):
	graph_data = filter_data(base_data, [(x_name, x_val)])
	grouped_data = graph_data.groupby(x_name)
	y = grouped_data.mean()[y_name].values
	val = 0.0
	if len(y) == 1:
		val = y[0]
	elif len(y) != 0:
		print("ERROR", y, len(y))
	return val


def metrics_real():
	file_append = 'real'
	filtered_data = data["metrics"].query("graph.str.contains('real')")
	for metric in metrics:
		fig, ax = plt.subplots()
		sns.barplot(x="graph",
					y=metric,
					hue="algo",
					palette="bright",
					data=filtered_data,
					ax=ax)
		sns.despine(left=True, ax=ax)
		ax.set(
			ylabel=metric_names[metric]["y_val"],
			xticklabels=[label.get_text()[5:] for label in ax.get_xticklabels()]
		)
		if metric is "time":
			ax.set_yscale('log')
		ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3)
		fig.suptitle(metric_names[metric]["description"] + " real graphs")
		plt.tight_layout(rect=(0, 0, 1, 0.96))
		fig.savefig(file_prefix + "metrics/" + metric_names[metric]["file_name"] + "_" + file_append + ".pdf")


def metrics_lfr():
	metrics_filter('LFR_om', 'om', algo_sets["clean"], "_clean")
	metrics_filter('LFR_om', 'om', algo_sets["ego"], "_ego")
	metrics_filter('LFR_mu', 'mu', algo_sets["ego"], "")


def set_xticklabels(ax, xlabel):
	if xlabel is "om":
		ax.set_xticklabels([str(i) for i in range(1, 6)])
	elif xlabel is "mu":
		ax.set_xticklabels([str(i) for i in range(0, 51, 5)])


def metrics_filter(graphs, xlabel, algos, name):
	filtered_data = data["metrics"].query("graph.str.contains(@graphs) and algo in @algos")
	if len(filtered_data) is 0:
		return
	for metric in metrics:
		fig, ax = plt.subplots()
		sns.lineplot(x="graph",
					 y=metric,
					 hue="algo",
					 ci=None,
					 style="algo",
					 markers=markers,
					 dashes=False,
					 linewidth=3,
					 markersize=8,
					 palette=colors,
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
		ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3)
		fig.suptitle(metric_names[metric]["description"] + ", LFR graphs" )
		plt.tight_layout(rect=(0, 0, 1, 0.96))
		fig.savefig(file_prefix + "metrics/" + metric_names[metric]["file_name"]
					+ "_" + graphs + name + ".pdf")


def num_comms_real():
	num_comms_filter('real')


def num_comms_lfr():
	num_comms_filter('LFR_om')
	num_comms_filter('LFR_mu')


def num_comms_filter(graphs):
	filtered_data = data["cover_num_comms"].query("graph.str.contains(@graphs)")
	fig, ax = plt.subplots()
	sns.lineplot(x="graph",
				 y="num_comms",
				 hue="algo",
				 style="algo",
				 ci=None,
				 # markers=True,
				 dashes=False,
				 linewidth=3,
				 # palette=sns.color_palette("hls", 13),
				 palette="bright",
				 data=filtered_data,
				 ax=ax)
	sns.despine(ax=ax)
	ax.set(ylim=(0, 300))
	ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
	fig.suptitle("Number of communities")
	plt.tight_layout(rect=(0, 0, 1, 0.96))
	fig.savefig(file_prefix + "communities/" + 'num_comms_' + graphs + '.pdf')
	# plt.close(fig)


# def comm_sizes():
# 	fig, ax = plt.subplots()
# 	graph = "Amazon_5000_no_small"
# 	markers = {
# 		"ground_truth": "x",
# 		"EgoSplitting_(PLM)": "o"
# 	}
# 	algo_names = data["comm_sizes"].groupby("algo").mean().index.values
#
# 	for algo_name in algo_names:
# 		filtered_data = filter_data(data["comm_sizes"], [("algo", algo_name), ("graph", graph)])
# 		grouped = filtered_data.groupby("comm_size").mean()
# 		x_val = grouped.index.values
# 		y_val = grouped["cnt"]
#
# 		ax.plot(x_val,
# 				y_val,
# 				# color=colors[algo_name],
# 				marker='x',#markers[algo_name],
# 				markersize=4,
# 				linestyle="None",
# 				label=algo_name,
# 				alpha=opacity
# 				)
#
# 	# Axes
# 	ax.set(xlabel="size",
# 		   ylabel="count",
# 		   title=graph + ", community sizes",
# 		   xlim=4,
# 		   ylim=0.5,
# 		   xscale='log',
# 		   yscale='log',
# 		   )
# 	ax.legend(loc="best")
#
# 	# Save plot
# 	fig.savefig(file_prefix + "comm_sizes.pdf")


def comm_sizes_real():
	comm_sizes_filter('real_')


def comm_sizes_lfr():
	comm_sizes_filter('LFR_om')
	comm_sizes_filter('LFR_mu')


def comm_sizes_filter(graphs):
	graph_names = sorted(list(set(data["cover_comm_sizes"].query("graph.str.contains(@graphs)")["graph"])))
	algos = algo_sets["comm_sizes"]

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
			palette=colors,
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


def comm_sizes_ego():
	comm_sizes_ego_filter('LFR_om')
	comm_sizes_ego_filter('LFR_mu')


def comm_sizes_ego_filter(graphs):
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


def node_comms_lfr():
	node_comms("LFR_om")
	node_comms("LFR_mu")


def node_comms(graphs):
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


def node_comms_ego():
	node_comms_ego_filter('LFR_om')
	node_comms_ego_filter('LFR_mu')


def node_comms_ego_filter(graphs):
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


def partition_counts():
	partition_counts_filter('LFR_om')
	partition_counts_filter('LFR_mu')


def partition_counts_filter(graphs):
	filtered_data = data["execution_info"].query("graph.str.contains(@graphs)"
												   " and info_name.str.contains('twoPlus')")
	fig, ax = plt.subplots()
	sns.lineplot(x="graph",
				 y="value",
				 hue="info_name",
				 style="algo",
				 markers=markers,
				 dashes=False,
				 linewidth=2,
				 markersize=8,
				 palette=colors,
				 data=filtered_data,
				 ax=ax)
	sns.despine(left=True, ax=ax)
	ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
	fig.suptitle("Partition counts per node (avg)")
	plt.tight_layout(rect=(0, 0, 1, 0.96))
	fig.savefig(file_prefix + 'partition_counts_' + graphs + '.pdf')
	# plt.close(fig)


def egonet_f1(graphs):
	filtered_data = data["execution_info"].query("graph.str.contains(@graphs)"
												 " and info_name.str.contains('egoF1Score')")
	fig, ax = plt.subplots()
	sns.lineplot(x="graph",
				 y="value",
				 hue="algo",
				 # style="algo",
				 ci=None,
				 markers=True,
				 dashes=False,
				 linewidth=2,
				 markersize=8,
				 palette="bright",
				 data=filtered_data,
				 ax=ax)
	sns.despine(left=True, ax=ax)
	ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
	fig.suptitle("Partition counts per node (avg)")
	plt.tight_layout(rect=(0, 0, 1, 0.96))
	fig.savefig(file_prefix + 'egonet_F1_' + graphs + '.pdf')
# plt.close(fig)


def egonet_comm_partition(graphs):
	filtered_data = data["ego_net_partitions"].query("graph.str.contains(@graphs)")

	# Nodes in other partitions (absolute)
	fig, ax = plt.subplots()
	sns.lineplot(
		x="comm_size",
		y="wrong_nodes",
		hue="algo",
		linewidth=3,
		palette=colors,
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
		palette=colors,
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
		palette=colors,
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


def egonet_partition_composition(graphs):
	filtered_data = data["ego_net_partition_composition"].query("graph.str.contains(@graphs)")

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
			palette=colors,
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
				hist_kws={"histtype": "step", "linewidth": 2,
				          "alpha": 1, "color": colors[algo_name]},
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


# metrics_real()
metrics_lfr()

# num_comms_real()
# num_comms_lfr()

# comm_sizes_lfr()
# comm_sizes_ego()

# node_comms_lfr()
# node_comms_ego()

# partition_counts()

# egonet_comm_partition("LFR_om_3")
# egonet_partition_composition("LFR_om_3")

# plt.show()