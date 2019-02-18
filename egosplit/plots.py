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
	data[name] = None
	filename = "cover_" + name + ".txt"
	try:
		file_data = pd.read_csv(filename, sep="\s+")
		data[name] = file_data
	except FileNotFoundError:
		print("File " + filename + " not found")

name = "partition_counts"
filename = name + '.txt'
try:
	data[name] = pd.read_csv(name + '.txt', sep="\s+")
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


# def plots_metrics(graph_names, file_append):
# 	graph_cnt = len(graph_names)
# 	algo_cnt = len(all_algo_names)
# 	bars_per_block = algo_cnt
# 	x_value_cnt = graph_cnt
#
# 	bar_width = 1 / (bars_per_block + 1.6)
# 	index = np.arange(x_value_cnt)
#
# 	for metric in metrics:
# 		fig, ax = plt.subplots()
# 		for i in range(algo_cnt):
# 			algo_name = all_algo_names[i]
# 			filtered_data = filter_data(data["metrics"], [("algo", algo_name)])
# 			y_val = []
# 			for graph_name in graph_names:
# 				y = get_value(filtered_data, "graph", graph_name, metric)
# 				y_val.append(y)
#
# 			ax.bar(index + i * bar_width,
# 				   y_val,
# 				   bar_width,
# 				   color=colors[algo_name],
# 				   label=algo_name,
# 				   error_kw=error_config,
# 				   alpha=opacity
# 				   )
#
# 		# Axes
# 		ax.set(ylabel=metric_names[metric]["y_val"],
# 			   title=metric_names[metric]["description"],
# 			   ylim=0,
# 			   xticks=index + bar_width * (bars_per_block - 1) / 2,
# 			   xticklabels=graph_names)
# 		ax.legend(loc="best")
#
# 		# Save plot
# 		fig.savefig(file_prefix + metric_names[metric]["file_name"] + "_" + file_append + ".pdf")


def metrics_real():
	file_append = 'real'
	filtered_data = data["metrics"].query("graph.str.contains('real')")
	for metric in metrics:
		sns.set(style="whitegrid", context="notebook", font_scale=0.8)
		fig, ax = plt.subplots()
		sns.barplot(x="graph",
					y=metric,
					hue="algo",
					data=filtered_data,
					ax=ax)
		sns.despine(left=True)
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
	metrics_filter('LFR_om')
	metrics_filter('LFR_mu')


def metrics_filter(graphs):
	filtered_data = data["metrics"].query("graph.str.contains(@graphs)")
	for metric in metrics:
		sns.set(style="whitegrid", context="notebook", font_scale=0.8)
		fig, ax = plt.subplots()
		sns.lineplot(x="graph",
					 y=metric,
					 hue="algo",
					 style="algo",
					 markers=True,
					 dashes=False,
					 # markersize=10,
					 linewidth=2,
					 markersize=8,
					 data=filtered_data,
					 ax=ax)
		# sns.despine(left=True)
		# print(ax.get_xticks())
		# print([label.get_text() for label in ax.get_xticklabels()])
		ax.set(
			ylabel=metric_names[metric]["y_val"],
			ylim=(-0.01, 1),
			# xticklabels=[label.get_text() for label in ax.get_xticklabels()]
		)
		ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3)
		fig.suptitle(metric_names[metric]["description"] + " LFR graphs" )
		plt.tight_layout(rect=(0, 0, 1, 0.96))
		fig.savefig(file_prefix + metric_names[metric]["file_name"]
					+ "_" + graphs + ".pdf")


def num_comms_real():
	filtered_data = data["num_comms"].query("graph.str.contains('real')")
	if len(filtered_data) == 0:
		return
	sns.set(style="whitegrid", context="notebook", font_scale=0.8)
	fig, ax = plt.subplots()
	sns.set(style="whitegrid", color_codes=True)
	sns.barplot(
		x="graph",
		y="num_comms",
		hue="algo",
		data=filtered_data,
		ax=ax,
		# xticklabels=[label.get_text()[5:] for label in ax.get_xticklabels()]
	)
	# sns.despine(left=True)
	ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3)
	fig.suptitle("Number of communities, real graphs")
	plt.tight_layout(rect=(0, 0, 1, 0.96))
	fig.savefig(file_prefix + 'num_comms.pdf')


def num_comms_lfr():
	num_comms_filter('LFR_om')
	num_comms_filter('LFR_mu')


def num_comms_filter(graphs):
	filtered_data = data["num_comms"].query("graph.str.contains(@graphs)")
	fig, ax = plt.subplots()
	sns.set(style="whitegrid")
	sns.lineplot(x="graph",
				 y="num_comms",
				 hue="algo",
				 style="algo",
				 markers=True,
				 dashes=False,
				 # markersize=10,
				 linewidth=2,
				 data=filtered_data,
				 ax=ax)
	# sns.despine(left=True)
	ax.set(title='Number of communities',
		   ylim=(0,200),
		   # yscale='log'
		   )
	fig.savefig(file_prefix + 'num_comms.pdf')


def comm_sizes():
	fig, ax = plt.subplots()
	graph = "Amazon_5000_no_small"
	markers = {
		"ground_truth": "x",
		"EgoSplitting_(PLM)": "o"
	}
	algo_names = data["comm_sizes"].groupby("algo").mean().index.values

	for algo_name in algo_names:
		filtered_data = filter_data(data["comm_sizes"], [("algo", algo_name), ("graph", graph)])
		grouped = filtered_data.groupby("comm_size").mean()
		x_val = grouped.index.values
		y_val = grouped["cnt"]

		ax.plot(x_val,
				y_val,
				# color=colors[algo_name],
				marker='x',#markers[algo_name],
				markersize=4,
				linestyle="None",
				label=algo_name,
				alpha=opacity
				)

	# Axes
	ax.set(xlabel="size",
		   ylabel="count",
		   title=graph + ", community sizes",
		   xlim=4,
		   ylim=0.5,
		   xscale='log',
		   yscale='log',
		   )
	ax.legend(loc="best")

	# Save plot
	fig.savefig(file_prefix + "comm_sizes.pdf")


def comm_sizes_sea(graphs):
	filtered_data = data["comm_sizes"].query("graph.str.contains(@graphs)")
	# for num_comms in range(5, 86, 10):
	# graph_name = "LFR_mu_" + str(num_comms).rjust(2, '0')
	fig, ax = plt.subplots()
	sns.set(style="whitegrid", color_codes=True)
	sns.violinplot(x="graph",
				   y="comm_size",
				   hue="algo",
				   scale="count",
				   data=filtered_data,
				   ax=ax)
	sns.despine(left=True)
	ax.set(title='Size of communities')
	fig.savefig(file_prefix + 'comm_sizes' + '.pdf')


def node_comms_sea():
	for num_comms in range(1, 6):
		graph_name = "LFR_GCE_" + str(num_comms)
		filtered_data = filter_data(data["node_comms"], [("graph", graph_name)])
		fig, ax = plt.subplots()
		sns.set(style="whitegrid", color_codes=True)
		sns.violinplot(x="graph",
					   y="num_comms",
					   hue="algo",
					   scale="count",
					   data=filtered_data,
					   ax=ax)
		sns.despine(left=True)
		ax.set(title='Number of communities per node')
		fig.savefig(file_prefix + 'node_comms' + str(num_comms) + '.pdf')


def comm_sizes_ego():
	comm_sizes_filter('LFR_om')
	comm_sizes_filter('LFR_mu')


def comm_sizes_filter(graphs):
	filtered_data = data["comm_sizes"].query(
		"(algo == 'EgoSplitting_(PLM)' or algo == 'ground_truth') and graph.str.contains(@graphs)")
	sns.set(style="whitegrid", context="notebook", font_scale=0.8)
	fig, ax = plt.subplots()
	sns.set(style="whitegrid", color_codes=True)
	sns.violinplot(x="graph",
				   y="comm_size",
				   hue="algo",
				   split=True,
				   scale="count",
				   data=filtered_data,
				   ax=ax)
	sns.despine(left=True)
	ax.set(
		ylabel="comm_size (log2)"
	)
	ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
	fig.suptitle("Size of coummunities, LFR graphs")
	plt.tight_layout(rect=(0, 0, 1, 0.96))
	fig.savefig(file_prefix + 'comm_sizes_' + graphs + '.pdf')


def node_comms_ego():
	node_comms_ego_filter('LFR_om')
	node_comms_ego_filter('LFR_mu')


def node_comms_ego_filter(graphs):
	filtered_data = data["node_comms"].query(
		"(algo == 'EgoSplitting_(PLM)' or algo == 'ground_truth') and graph.str.contains(@graphs)")
	sns.set(style="whitegrid", context="notebook", font_scale=0.8)
	fig, ax = plt.subplots()
	sns.set(style="whitegrid", color_codes=True)
	sns.violinplot(x="graph",
				   y="num_comms",
				   hue="algo",
				   split=True,
				   scale="count",
				   data=filtered_data,
				   ax=ax)
	sns.despine(left=True)
	ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3)
	fig.suptitle("Number of communities per node, LFR graphs")
	plt.tight_layout(rect=(0, 0, 1, 0.96))
	fig.savefig(file_prefix + 'node_comms_' + graphs + '.pdf')


def partition_counts():
	# partition_counts_filter('LFR_om')
	partition_counts_filter('LFR_mu')


def partition_counts_filter(graphs):
	filtered_data = data["partition_counts"].query("graph.str.contains(@graphs)")
	sns.set(style="whitegrid", context="notebook", font_scale=0.8)
	fig, ax = plt.subplots()
	sns.lineplot(x="graph",
				 y="count",
				 hue="partition_name",
				 style="algo",
				 markers=True,
				 dashes=False,
				 linewidth=2,
				 markersize=8,
				 data=filtered_data,
				 ax=ax)
	sns.despine(left=True)
	ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, prop={'size': 9})
	fig.suptitle("Partition counts per node (avg)")
	plt.tight_layout(rect=(0, 0, 1, 0.96))
	fig.savefig(file_prefix + 'partition_counts' + '.pdf')


# comm_sizes_sea('LFR_mu')
# node_comms_sea()

comm_sizes_ego()
node_comms_ego()

# metrics_real()
# metrics_lfr()

num_comms_real()
num_comms_lfr()

partition_counts()

# plt.show()