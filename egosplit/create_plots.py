from plot_scripts.comm_sizes import create_comm_sizes_plots
from plot_scripts.comms_per_node import create_node_comms_plots
from plot_scripts.egonet import create_egonet_plots
from plot_scripts.read_data import read_data
from plot_scripts.config import set_sns_style
from plot_scripts.base_plot import make_plot
from plot_scripts.comm_sizes import metric_names


print("Reading results...")
data = read_data()
set_sns_style()

print("Creating plots...")
plots = [
	"metrics",
	"comm_sizes",
	# "num_comms",
	# "ego_net_partitioning",
]


# *****************************************************************************
# *                                  Metrics                                  *
# *****************************************************************************
if "metrics" in plots:
	metrics = [
		'time',
		'nmi',
		# 'f1',
		# 'f1_rev',
		# 'entropy',
		# 'entropy2',
		# 'entropy3',
		# 'entropy4',
	]
	for metric in metrics:
		# break
		make_plot(
			data=data["metrics"],
			graphs="LFR_om",
			xlabel="om",
			algo_match="",
			# add_algos=["OSLOM"],
			remove_algo_part="Ego_PLM_Info_",
			title=metric_names[metric]["description"],
			file_name="metrics/" + metric_names[metric]["file_name"],
			x="graph",
			y=metric,
			hue="algo",
			plot_args={
				"dashes": False,
				# "markers": False,
			},
			ax_set={
				"ylim": metric_names[metric]["ylim"],
				"ylabel": metric_names[metric]["y_val"]
			}
		)
		# make_plot(
		# 	data=data["metrics"],
		# 	graphs="LFR_om",
		# 	xlabel="om",
		# 	algo_match="_triangles",
		# 	title=metric_names[metric]["description"] + ", extend(triangles)",
		# 	file_name="metrics/" + metric_names[metric]["file_name"] + "_triangles",
		# 	x="graph",
		# 	y=metric,
		# 	hue="algo",
		# 	plot_args={
		# 	},
		# 	ax_set={
		# 		"ylim": metric_names[metric]["ylim"],
		# 		"ylabel": metric_names[metric]["y_val"]
		# 	}
		# )

	# metric = "nmi"
	# x = {"name": "factor", "create_from": "algo", "str_start": "_f-", "str_end": "*"}
	# hue = {"name": "exponent", "create_from": "algo", "str_start": "*e-", "str_end": ""}
	# make_plot(
	# 	data=data["metrics"],
	# 	graphs="LFR_om",
	# 	algo_match="Ego_PLM_",
	# 	title=metric_names[metric]["description"],
	# 	file_name="2_dim_metrics/" + metric_names[metric]["file_name"],
	# 	x=x,
	# 	y=metric,
	# 	hue=hue,
	# 	plot_args={
	#
	# 	},
	# 	ax_set={
	# 		"ylim": metric_names[metric]["ylim"],
	# 		"ylabel": metric_names[metric]["y_val"]
	# 	}
	# )


# *****************************************************************************
# *                                  Comm sizes                               *
# *****************************************************************************
if "comm_sizes" in plots:
	make_plot(
		plot_type="swarm",
		data=data["cover_comm_sizes"],
		graphs="LFR_om_",
		xlabel="",
		algo_match="",
		# add_algos=["Ground_Truth"],
		remove_algo_part=["Ego_PLM_Info_"],
		title="Community Sizes",
		file_name="communities/" + "comm_sizes",
		x="algo",
		hue="algo",
		y="comm_size",
		plot_args={
		},
		ax_set={
			"ylim": (0, 10),
			"ylabel": "size (log2)",
			"xticklabels": [""],
		}
	)


# *****************************************************************************
# *                            Number of communities                          *
# *****************************************************************************
if "num_comms" in plots:
	make_plot(
		data=data["cover_num_comms"],
		graphs="LFR_om",
		xlabel="om",
		algo_match="Ego_",
		add_algos=["Ground_Truth"],
		remove_algo_part="Ego_PLM_",
		title="Number of communities" + ", PLM(1.0)",
		file_name="communities/" + "num_comms",
		x="graph",
		y="num_comms",
		hue="algo",
		plot_args={
			"dashes": False,
			"markers": False,
		},
		ax_set={
			"ylim": 0,
		}
	)


# *****************************************************************************
# *                            Ego-Net partitioning                           *
# *****************************************************************************
if "ego_net_partitioning" in plots:
	ego_metrics = [
		"community_cohesion",
		"partition_exclusivity",
		"ego_partitioning_score",
		"merged_external_nodes",
		"parts_per_comm",
		"comms_per_part",
	]
	algos = [
		# {"filter": "_edges", "title": "extend (edges)", "file": "_edges"},
		# {"filter": "_triangles", "title": "extend (triangles)", "file": "_triangles"},
		# {"filter": "Ego_PLM_1.0_", "title": "PLM(1.0)", "file": "_PLM_1.0"},
		{"filter": "", "title": "", "file": ""},
	]
	for ego_metric in ego_metrics:
		# break
		for algo in algos:
			make_plot(
				data=data["ego_net_metrics"].query("metric_name in @ego_metric"),
				graphs="LFR_om",
				xlabel="om",
				algo_match=algo["filter"],
				remove_algo_part="Ego_",
				title="Ego-Net Metrics, " + ego_metric + ", " + algo["title"],
				file_name="ego_partition/metrics/" + ego_metric + algo["file"],
				x="graph",
				y="value",
				hue="algo",
				plot_args={
					# "style": "metric_name",
					"dashes": False,
					# "markers": False,
				},
				ax_set={
					"ylim": (0, 1.05),
				}
			)

	# Algo parameters on x-axis
	for ego_metric in ego_metrics:
		break
		x = {"name": "factor", "create_from": "algo", "str_start": "_f-", "str_end": "*"}
		hue = {"name": "exponent", "create_from": "algo", "str_start": "*e-", "str_end": ""}
		make_plot(
			data=data["ego_net_metrics"].query("metric_name in @ego_metric"),
			graphs="LFR_om",
			algo_match="Ego_PLM_",
			title="Ego-Net Metrics, " + ego_metric,
			file_name="ego_partition/2_dim_metrics/" + ego_metric,
			x=x,
			y="value",
			hue=hue,
			plot_args={
			},
			ax_set={
				"ylim": (0, 1.05),
			}
		)

	# Metrics per Ego-Net
	for ego_metric in ego_metrics:
		# break
		for algo in algos:
			make_plot(
				data=data["ego_net_ego_metrics"].query("metric_name in @ego_metric"),
				graphs="LFR_om",
				# xlabel="om",
				algo_match=algo["filter"],
				title="Ego-Net Metrics, " + ego_metric + ", " + algo["title"],
				file_name="ego_partition/ego_metrics/" + ego_metric + algo["file"],
				x="ego_net_size",
				y="value",
				hue="algo",
				plot_args={
					# "style": "metric_name",
					"dashes": False,
					# "markers": False,
					"markersize": 3,
				},
				ax_set={
					"ylim": (0, 1.05),
				}
			)


# # Create plots
# create_metric_plots(data)
#
# create_num_comms_plots(data)
#
# create_comm_sizes_plots(data)
#
# create_node_comms_plots(data)
#
# create_partition_counts_plots(data)
#
# create_egonet_plots(data)
# plt.show()
