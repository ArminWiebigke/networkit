import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import *


def metrics_lfr(data):
	# metrics_filter('LFR_om', 'om', algo_sets["clean"], "_clean")
	# metrics_filter('LFR_om', 'om', algo_sets["ego"], "_ego")
	# metrics_filter('LFR_mu', 'mu', algo_sets["ego"], "")
	# metrics_filter('LFR_om', 'om', algo_sets["ego_parameters"], "_parameters")
	metrics_filter(data, 'LFR_om', 'om', "edgeScore_0.4", "_0.4", True)
	metrics_filter(data, 'LFR_om', 'om', "edgeScore_0.6", "_0.6", True)
	metrics_filter(data, 'LFR_om', 'om', "edgeScore_0.8", "_0.8", True)
	metrics_filter(data, 'LFR_om', 'om', "edgeScore_1.0", "_1.0", True)


def metrics_filter(data, graphs, xlabel, algos, name, match_algo_name=False):
	filtered_data = data["metrics"].query("graph.str.contains(@graphs)")
	if match_algo_name:
		filtered_data = filtered_data.query("algo.str.contains(@algos)")
	else:
		filtered_data = filtered_data.query("algo in @algos")
	if len(filtered_data) is 0:
		return

	for metric in ["nmi"]:#metrics:
		fig, ax = plt.subplots()
		sns.lineplot(x="graph",
		             y=metric,
		             hue="algo",
		             ci=None,
		             style="algo",
		             # markers=True,
		             # markers=markers,
		             dashes=False,
		             linewidth=2,
		             markersize=6,
		             palette=sns.light_palette("green", 20),
		             # palette=colors,
		             data=filtered_data,
		             ax=ax)
		sns.despine(ax=ax)
		ax.set(
			xlabel=xlabel,
			ylabel=metric_names[metric]["y_val"],
		)
		legend_handles, legend_labels = ax.get_legend_handles_labels()
		legend_labels = [l[18:] for l in legend_labels]

		set_xticklabels(ax, xlabel)
		if metric != "time":
			ax.set(ylim=(-0.01, 1))
		set_layout(ax, legend_handles, legend_labels)
		fig.suptitle(metric_names[metric]["description"] + ", LFR graphs" )
		file_name = file_prefix + "metrics/" + metric_names[metric]["file_name"] + "_" \
		            + graphs + name
		fig.savefig(file_name + ".pdf")
