import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# Read input *******************************************************************
data = dict()
data["metrics"] = None
filename = "results.txt"
try:
	file_data = pd.read_csv(filename, sep="\s+")
	data["metrics"] = file_data
except FileNotFoundError:
	print("File " + filename + " not found")

metrics = ["F1", "F1_rev", "NMI", "time"]
all_graph_names = data["metrics"].groupby("graph").mean().index.values
all_algo_names = data["metrics"].groupby("algo").mean().index.values

for name in ["comm_sizes", "node_comms", "num_comms"]:
	filename = "cover_" + name + ".txt"
	try:
		data[name] = pd.read_csv(filename, sep="\s+")
	except FileNotFoundError:
		print("File " + filename + " not found")

result_dir = "results/"
for name in ["execution_info", "ego_net_partitions", "ego_net_partition_composition"]:
	filename = name + '.result'
	try:
		data[name] = pd.read_csv(result_dir + filename, sep="\s+")
	except FileNotFoundError:
		print("File " + filename + " not found")


# Create Plots *****************************************************************
file_prefix = "plots/"
opacity = 0.7
error_config = {'ecolor': '0.3',
				'capsize': 4}
colors = {
	"ground_truth": "xkcd:black",
	"EgoSplitting_(PLP)": "xkcd:goldenrod",
	"EgoSplitting_(PLM)": "xkcd:red",
	"EgoSplitting_(LPPotts)": "xkcd:light orange",
	"EgoSplitting_(LPPotts_par)": "xkcd:orange",
	"OLP": "xkcd:green",
	"GCE": "xkcd:blue",
	"MOSES": "xkcd:dark teal",
	"OSLOM": "xkcd:fuchsia",
}
metric_names = {
	"F1": {
		"description": "F1 Score",
		"y_val": "F1",
		"file_name": "F1_score"
	},
	"F1_rev": {
		"description": "F1 Score (reversed)",
		"y_val": "F1",
		"file_name": "F1_score_rev"
	},
	"NMI": {
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
		fig.savefig(file_prefix + metric_names[metric]["file_name"] + "_" + file_append + ".pdf")


def metrics_lfr():
	metrics_filter('LFR_om', 'om')
	metrics_filter('LFR_mu', 'mu')


def metrics_filter(graphs, xlabel):
	filtered_data = data["metrics"].query("graph.str.contains(@graphs)")
	for metric in metrics:
		fig, ax = plt.subplots()
		sns.lineplot(x="graph",
					 y=metric,
					 hue="algo",
					 ci=None,
					 style="algo",
					 markers=True,
					 dashes=False,
					 linewidth=3,
					 markersize=8,
					 palette="bright",
					 data=filtered_data,
					 ax=ax)
		sns.despine(ax=ax)
		ax.set(
			xlabel=xlabel,
			ylabel=metric_names[metric]["y_val"],
			# xticklabels=[label.get_text() for label in ax.get_xticklabels()]
		)
		if metric != "time":
			ax.set(ylim=(-0.01, 1))
		ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3)
		fig.suptitle(metric_names[metric]["description"] + " LFR graphs" )
		plt.tight_layout(rect=(0, 0, 1, 0.96))
		fig.savefig(file_prefix + metric_names[metric]["file_name"]
					+ "_" + graphs + ".pdf")


def num_comms_real():
	num_comms_filter('real')


def num_comms_lfr():
	num_comms_filter('LFR_om')
	num_comms_filter('LFR_mu')


def num_comms_filter(graphs):
	filtered_data = data["num_comms"].query("graph.str.contains(@graphs)")
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
	fig.savefig(file_prefix + 'num_comms_' + graphs + '.pdf')
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
	graph_names = sorted(list(set(data["comm_sizes"].query("graph.str.contains(@graphs)")["graph"])))
	for graph_name in graph_names:
		filtered_data = data["comm_sizes"].query(("graph == @graph_name"))
		fig, ax = plt.subplots()
		sns.violinplot(x="graph",
					   y="comm_size",
					   hue="algo",
					   scale="count",
					   linewidth=0.5,
					   inner=None,
					   palette="bright",
					   data=filtered_data,
					   ax=ax,
					   )
		sns.despine(ax=ax)
		ax.set(ylabel="size (log2)",
			   ylim=(0,10))
		ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
		fig.suptitle("Size of communities")
		plt.tight_layout(rect=(0, 0, 1, 0.96))
		fig.savefig(file_prefix + 'comm_sizes_' + graph_name + '.pdf')
		# plt.close(fig)


def comm_sizes_ego():
	comm_sizes_ego_filter('LFR_om')
	comm_sizes_ego_filter('LFR_mu')


def comm_sizes_ego_filter(graphs):
	filtered_data = data["comm_sizes"].query(
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
	fig.savefig(file_prefix + 'comm_sizes_' + graphs + '_ego.pdf')
	# plt.close(fig)


def node_comms_lfr():
	node_comms("LFR_om")
	node_comms("LFR_mu")


def node_comms(graphs):
	graph_names = sorted(list(set(data["node_comms"].query("graph.str.contains(@graphs)")["graph"])))
	for graph_name in graph_names:
		filtered_data = data["node_comms"].query(("graph == @graph_name"))
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
		fig.savefig(file_prefix + 'node_comms_' + graph_name + '.pdf')
		# plt.close(fig)


def node_comms_ego():
	node_comms_ego_filter('LFR_om')
	node_comms_ego_filter('LFR_mu')


def node_comms_ego_filter(graphs):
	filtered_data = data["node_comms"].query(
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
	fig.savefig(file_prefix + 'node_comms_' + graphs + '.pdf')
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
	sns.lineplot(
		x="comm_size",
		y="wrong_nodes",
		hue="algo",
		data=filtered_data,
	)

	# sns.catplot(
	# 	x="comm_size",
	# 	y="wrong_nodes",
	# 	kind="swarm",
	# 	hue="num_partitions",
	# 	col="algo",
	# 	data=filtered_data,
	# )

	for algo in filtered_data.groupby("algo").mean().index.values:
		sns.jointplot(
			x="comm_size",
			y="wrong_nodes",
			kind="hex",
			data=filtered_data.query("algo == @algo"),
			xlim=(0, 38),
			ylim=(-1, 20),
		)
		plt.suptitle(algo)
		plt.tight_layout(rect=(0, 0, 1, 0.96))

	for algo in filtered_data.groupby("algo").mean().index.values:
		sns.jointplot(
			x="comm_size",
			y="wrong_nodes",
			kind="kde",
			data=filtered_data.query("algo == @algo"),
			xlim=(0, 38),
			ylim=(-1, 20),
		)
		plt.suptitle(algo)
		plt.tight_layout(rect=(0, 0, 1, 0.96))


def egonet_partition_composition(graphs):
	filtered_data = data["ego_net_partition_composition"].query("graph.str.contains(@graphs)")
	# fig, ax = plt.subplots()
	# sns.lineplot(
	# 	x="partition_size",
	# 	y="wrong_nodes",
	# 	hue="algo",
	# 	ci=None,
	# 	data=filtered_data,
	# 	ax=ax,
	# )
	# fig, ax = plt.subplots()
	# sns.lineplot(
	# 	x="partition_size",
	# 	y="wrong_percentage",
	# 	hue="algo",
	# 	ci=None,
	# 	data=filtered_data,
	# 	ax=ax,
	# )

	fig, ax = plt.subplots()
	sns.violinplot(x="graph",
	               y="partition_size",
	               hue="algo",
	               scale="count",
	               linewidth=1,
	               inner=None,
	               palette="bright",
	               data=filtered_data,
	               ax=ax)

	# for algo in filtered_data.groupby("algo").mean().index.values:
	# 	sns.jointplot(
	# 		x="partition_size",
	# 		y="wrong_nodes",
	# 		kind="hex",
	# 		data=filtered_data.query("algo == @algo"),
	# 		# xlim=(0, 80),
	# 		# ylim=(-1, 50),
	# 	)
	# 	plt.suptitle(algo)
	# 	plt.tight_layout(rect=(0, 0, 1, 0.96))
	#
	# for algo in filtered_data.groupby("algo").mean().index.values:
	# 	sns.jointplot(
	# 		x="partition_size",
	# 		y="wrong_nodes",
	# 		kind="kde",
	# 		data=filtered_data.query("algo == @algo"),
	# 		# xlim=(0, 80),
	# 		# ylim=(-1, 50),
	# 	)
	# 	plt.suptitle(algo)
	# 	plt.tight_layout(rect=(0, 0, 1, 0.96))


# metrics_real()
# metrics_lfr()

# num_comms_real()
# num_comms_lfr()

comm_sizes_lfr()
# comm_sizes_ego()

node_comms_lfr()
# node_comms_ego()

# partition_counts()

# egonet_f1('LFR_om')
# egonet_f1('LFR_mu')

# egonet_comm_partition("LFR_om_3")
# egonet_partition_composition("LFR_om_3")

plt.show()