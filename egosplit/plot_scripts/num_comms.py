import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import *


def num_comms_lfr(data):
	num_comms_filter(data, 'LFR_om')


def num_comms_filter(data, graphs):
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
	fig.suptitle("Number of communities")
	set_layout(ax)
	fig.savefig(file_prefix + "communities/" + 'num_comms_' + graphs + '.pdf')
# plt.close(fig)
