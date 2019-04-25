import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from collections import defaultdict

from .config import *
from .utils import filter_data


def create_num_comms_plots(data):
	num_comms_filter(data, 'LFR_om', "om", "", 'Ego_PLM_1.0_', "")
	pass


def num_comms_filter(data, graphs, xlabel, algos, remove_algo_part, name):
	filtered_data = filter_data(data["cover_num_comms"], graphs, algos)
	if len(filtered_data) is 0:
		return
	fig, ax = plt.subplots()
	sns.lineplot(x="graph",
	             y="num_comms",
	             hue="algo",
	             style="algo",
	             ci=None,
	             markers=defaultdict(lambda: "x"),
	             dashes=False,
	             linewidth=3,
	             # palette=sns.color_palette("hls", 13),
	             palette="bright",
	             data=filtered_data,
	             ax=ax)
	sns.despine(ax=ax)
	# ax.set(ylim=(0, 300))
	fig.suptitle("Number of communities")
	legend_handles, legend_labels = ax.get_legend_handles_labels()
	if type(algos) == str and remove_algo_part == "":
		remove_algo_part = algos
	legend_labels = [l.replace(remove_algo_part, '') for l in legend_labels]
	set_layout(ax, legend_handles, legend_labels)
	fig.savefig(file_prefix + "communities/" + 'num_comms_' + graphs + '.pdf')
	# plt.close(fig)
