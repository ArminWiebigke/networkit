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


def set_legend(ax, num_columns, legend_handles, legend_labels):
	legend_args = {
		"fontsize": legend_font_size(),
		"ncol": num_columns,
		"loc": "lower center",
		"bbox_to_anchor": (0.5, 1.01),
	}
	ax.legend(**legend_args, handles=legend_handles, labels=legend_labels)


def set_layout():
	# plt.tight_layout(rect=(0, 0, 1, 0.96))
	plt.tight_layout(pad=0.5)


def legend_font_size():
	return 10
