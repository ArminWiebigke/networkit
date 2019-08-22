from collections import defaultdict

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import pyplot as plt

import egosplit.benchmarks.evaluation.benchmark_metric as bm

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
	bm.F1.get_name(): {
		"ylabel": "F1 Score",
		"y_val": "F1-Score",
		"file_name": "F1",
		"ylim": (0, 1.05),
	},
	bm.F1_rev.get_name(): {
		"ylabel": "F1 Score (reversed)",
		"y_val": "F1-Score (reversed)",
		"file_name": "F1_rev",
		"ylim": (0, 1.05),
	},
	bm.NMI.get_name(): {
		"ylabel": "NMI",
		"y_val": "NMI",
		"file_name": "NMI",
		"ylim": (0, 1.05),
	},
	bm.Time.get_name(): {
		"ylabel": r"Running Time / (n + m) [\SI{}{\micro\second}]",
		"y_val": "time / (n + m)",
		"file_name": "time",
		"ylim": 0,
	},
}

ego_metric_ylim = defaultdict(lambda: 0)
ego_metric_ylim['Persona Recall'] = (0, 1.05)

algo_sets = dict()
algo_sets["clean"] = [
	"Ego_PLM", "Ego_PLM_clean_OSLOM", "OSLOM",
]
algo_sets["ego"] = [
	"Ego_PLP", "Ego_LPPotts", "Ego_PLM", "Ego_Infomap", "Ego_Surprise",
]
algo_sets["comm_sizes"] = [
	"ground_truth", "Ego_PLM", "Ego_Infomap", "Ego_Surprise", "Ego_PLP", "MOSES",
	"OSLOM", "GCE"
]
algo_sets["ego"] = [
	"Ego_PLM", "Ego_Infomap", "Ego_PLP", "Ego_Surprise"
]
algo_sets["ego_parameters"] = [
	"Ego_PLM_base",
	"Ego_PLM_base_add",
	"Ego_PLM_simpleNN",
	"Ego_PLM_simpleNN_add",
	"Ego_PLM_edgeScores",
	"Ego_PLM_edgeScores_add",
]


def set_sns_style():
	plt.rc('text', usetex=True)
	preamble = \
		r'\usepackage{siunitx} ' \
		r'\usepackage[utf8]{inputenc}' \
		r'\usepackage[T1]{fontenc}' \
		r'\usepackage{upgreek}' \
		r'\usepackage{lmodern}'
	plt.rc('text.latex', preamble=preamble)
	sns.set(context="notebook", style="whitegrid", palette="bright",
	        font_scale=1.8,
	        )


def set_legend(ax, legend_handles=None, legend_labels=None):
	legend_args = {
		**get_legend_args(),
		"loc": "lower center",
		"bbox_to_anchor": (0.5, 1.01)
	}
	if legend_handles is None:
		ax.legend(**legend_args)
	else:
		ax.legend(**legend_args, handles=legend_handles, labels=legend_labels)


def get_legend_args():
	legend_args = {
		# "ncol": 3,
		"prop": {'size': 8}
	}
	return legend_args


def set_layout():
	plt.tight_layout(rect=(0, 0, 1, 0.96))
