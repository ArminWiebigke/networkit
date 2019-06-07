import matplotlib.pyplot as plt
import seaborn as sns

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
		"file_name": "F1_score",
		"ylim": (0, 1.05),
	},
	"f1_rev": {
		"description": "F1 Score (reversed)",
		"y_val": "F1",
		"file_name": "F1_score_rev",
		"ylim": (0, 1.05),
	},
	"nmi": {
		"description": "NMI Score",
		"y_val": "NMI",
		"file_name": "NMI_score",
		"ylim": (0, 1.05),
	},
	"time": {
		"description": "Running Time",
		"y_val": "time (s)",
		"file_name": "time",
		"ylim": 0,
	},
	"entropy": {
		"description": "Entropy",
		"y_val": "Entropy",
		"file_name": "entropy",
		"ylim": 0,
	},
	"entropy2": {
		"description": "Entropy2",
		"y_val": "Entropy",
		"file_name": "entropy2",
		"ylim": 0,
	},
	"entropy3": {
		"description": "Entropy3",
		"y_val": "Entropy",
		"file_name": "entropy3",
		"ylim": 0,
	},
	"entropy4": {
		"description": "Entropy4",
		"y_val": "Entropy",
		"file_name": "entropy4",
		"ylim": 0,
	},
}

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
	sns.set(context="notebook", style="whitegrid", palette="bright", font_scale=0.8)


def set_xticklabels(ax, xlabel):
	labels = ax.get_xticklabels()
	def cut_label(label, xlabel):
		begin = label.find(xlabel) + len(xlabel)
		end = label.find('_', begin)
		if end == -1:
			return label[begin:]
		return label[begin:end]

	labels = [cut_label(l.get_text(), xlabel) for l in labels]
	ax.set_xticklabels(labels)


def set_layout(ax, legend_handles=None, legend_labels=None):
	legend_args = {
		"loc": "lower center",
		"bbox_to_anchor": (0.5, 1.01),
		"ncol": 3,
		"prop": {'size': 9}
	}
	if legend_handles is None:
		ax.legend(**legend_args)
	else:
		ax.legend(**legend_args, handles=legend_handles, labels=legend_labels)
	plt.tight_layout(rect=(0, 0, 1, 0.96))