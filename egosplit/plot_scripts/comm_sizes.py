import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import *
from .utils import filter_data


def create_comm_sizes_plots(data):
	pass



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
