import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Read input *******************************************************************
data = None
filename = "results.txt"
try:
	file_data = pd.read_csv(filename, sep="\s+")
	data = file_data
except FileNotFoundError:
	print("File " + filename + " not found")

metrics = ["F1", "F1_(rev)", "NMI"]
all_graph_names = data.groupby("graph").mean().index.values
all_algo_names = data.groupby("algo").mean().index.values

# Create Plots *****************************************************************
file_prefix = "plots/"
opacity = 0.7
error_config = {'ecolor': '0.3',
				'capsize': 4}
colors = {
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
	"F1_(rev)": {
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


def plots_large():
	graph_names = [g for g in all_graph_names if not "LFR" in g]
	plots_metrics(graph_names, "large")
	plots_time(graph_names, "large")


def plots_lfr():
	graph_names = [g for g in all_graph_names if "LFR" in g]
	plots_metrics(graph_names, "lfr")
	plots_time(graph_names, "lfr")


def plots_metrics(graph_names, file_append):
	graph_cnt = len(graph_names)
	algo_cnt = len(all_algo_names)
	bars_per_block = algo_cnt
	x_value_cnt = graph_cnt

	bar_width = 1 / (bars_per_block + 1.6)
	index = np.arange(x_value_cnt)

	for metric in metrics:
		fig, ax = plt.subplots()
		for i in range(algo_cnt):
			algo_name = all_algo_names[i]
			filtered_data = filter_data(data, [("algo", algo_name)])
			y_val = []
			for graph_name in graph_names:
				y = get_value(filtered_data, "graph", graph_name, metric)
				y_val.append(y)

			ax.bar(index + i * bar_width,
				   y_val,
				   bar_width,
				   color=colors[algo_name],
				   label=algo_name,
				   error_kw=error_config,
				   alpha=opacity
				   )

		# Axes
		ax.set(ylabel=metric_names[metric]["y_val"],
			   title=metric_names[metric]["description"],
			   ylim=0,
			   xticks=index + bar_width * (bars_per_block - 1) / 2,
			   xticklabels=graph_names)
		ax.legend(loc="best")

		# Save plot
		fig.savefig(file_prefix + metric_names[metric]["file_name"] + "_" + file_append + ".pdf")


def plots_time(graph_names, file_append):
	graph_cnt = len(graph_names)
	algo_cnt = len(all_algo_names)
	bars_per_block = algo_cnt
	x_value_cnt = graph_cnt

	bar_width = 1 / (bars_per_block + 1.6)
	index = np.arange(x_value_cnt)

	metric = "time"
	fig, ax = plt.subplots()
	for i in range(algo_cnt):
		algo_name = all_algo_names[i]
		filtered_data = filter_data(data, [("algo", algo_name)])
		y_val = []
		for graph_name in graph_names:
			y = get_value(filtered_data, "graph", graph_name, metric)
			y_val.append(y)

		ax.bar(index + i * bar_width,
			   y_val,
			   bar_width,
			   color=colors[algo_name],
			   label=algo_name,
			   error_kw=error_config,
			   alpha=opacity
			   )

	# Axes
	ax.set(ylabel=metric_names[metric]["y_val"],
		   title=metric_names[metric]["description"],
		   xticks=index + bar_width * (bars_per_block - 1) / 2,
		   xticklabels=graph_names,
		   yscale='log',
		   )
	ax.legend(loc="best")

	# Save plot
	fig.savefig(file_prefix +metric_names[metric]["file_name"] + "_" + file_append + ".pdf")


plots_large()
plots_lfr()
# plt.show()