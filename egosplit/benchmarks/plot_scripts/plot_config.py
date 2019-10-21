import seaborn as sns
from matplotlib import pyplot as plt


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
	# plt.tight_layout(rect=(0, 0, 1, 0.96))
	plt.tight_layout(pad=0.5)
