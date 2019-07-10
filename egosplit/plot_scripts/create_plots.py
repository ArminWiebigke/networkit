from matplotlib.ticker import MultipleLocator

from .config import metric_names
from .draw_plot import make_plot, PlotType


algo_matches = ["Ego", "GCE"]
remove_algo_part = ["Ego_", "eiden_Mod", "_e!"]


def run(data):
	graph_sets = {
		"om": {
			'graph_filter': "_om_",
			'x': {"name": "communities per node", "create_from": "graph", "str_start": "_om_",
			      "str_end": "_", 'type': float},
			'ax_set': {
				# 'xticks': [1, 1.2, 1.4, 1.6, 1.8, 2, 3, 4, 5],
				# 'minor_locator': (MultipleLocator(0.2)),
			}
		},
		"mu": {
			'graph_filter': "_mu_",
			'x': {"name": "mixing factor", "create_from": "graph", "str_start": "_mu_",
			      "str_end": "_", 'type': float},
			'ax_set': {
				# 'xticklabels': [1, 1.2, 1.4, 1.6, 1.8, 2, 3, 4, 5],
			}
		}
	}
	for graph_set_name, graph_set_params in graph_sets.items():
		params = (data, graph_set_name, graph_set_params)
		metric_plots(*params)
		comm_sizes_plots(*params)
		num_comms_plots(*params)
		ego_net_plots(*params)
		ego_net_plots_per_graph(*params)


# *****************************************************************************
# *                                  Metrics                                  *
# *****************************************************************************
def metric_plots(data, graph_set_name, graph_set_params):
	print("Plots for Metrics")
	metrics = [
		'time',
		'nmi',
		'f1',
		'f1_rev',
		# 'entropy',
		# 'entropy2',
		# 'entropy3',
		# 'entropy4',
	]
	for metric in metrics:
		# break
		make_plot(
			# plot_type=PlotType.bar,
			data=data["metrics"],
			graph_filter=graph_set_params["graph_filter"],
			# xlabel="om",
			algo_matches=algo_matches,
			# add_algos=["GCE_1.0", "GCE_1.5"],
			remove_algo_part=remove_algo_part,
			title=metric_names[metric]["description"],
			file_name="metrics/{}_{}".format(metric_names[metric]["file_name"],
			                                 graph_set_name),
			one_plot_per_graph=False,
			x=graph_set_params['x'],
			y=metric,
			hue="algo",
			plot_args={
				"ci": "sd",
			},
			ax_set={
				**graph_set_params['ax_set'],
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
def comm_sizes_plots(data, graph_set_name, graph_set_params):
	print("Plots for comm sizes")
	make_plot(
		plot_type=PlotType.swarm,
		data=data["cover_comm_sizes"],
		graph_filter=graph_set_params["graph_filter"],
		algo_matches=algo_matches,
		add_algos=["Ground_Truth"],
		remove_algo_part=remove_algo_part,
		title="Community Sizes",
		file_name="communities/comm_sizes_{}".format(graph_set_name),
		one_plot_per_graph=True,
		x=graph_set_params['x'],
		hue="algo",
		y="comm_size",
		plot_args={
			# "size": 2.5,
			"dodge": True,
		},
		ax_set={
			# "ylim": (2, 8),
			# "ylim": 2,
			"ylabel": "size (log2)",
			# "xticklabels": [""],
		}
	)


# *****************************************************************************
# *                            Number of communities                          *
# *****************************************************************************
def num_comms_plots(data, graph_set_name, graph_set_params):
	print("Plots for num comms")
	make_plot(
		# plot_type=PlotType.bar,
		data=data["cover_num_comms"],
		graph_filter=graph_set_params["graph_filter"],
		algo_matches=algo_matches,
		add_algos=["Ground_Truth"],
		remove_algo_part=remove_algo_part,
		title="Number of communities",
		file_name="communities/num_comms_{}".format(graph_set_name),
		one_plot_per_graph=False,
		x=graph_set_params['x'],
		y="num_comms",
		hue="algo",
		plot_args={
		},
		ax_set={
			"ylim": 0,
		}
	)


# *****************************************************************************
# *                             Ego-Net partition                             *
# *****************************************************************************
ego_metrics = [
	# "community_cohesion",
	# "partition_exclusivity",
	# "ego_partition_score",
	# "merged_external_nodes",
	# "parts_per_comm",
	# "comms_per_part",
	"extended_nodes",
	"external_nodes",
	"external_nodes_added",
	"external_nodes_added_total",
	"conductance",
	"conductance_ratio",
	"intra_edges",
	"intra_ratio",
	"num_components",
	"separate_nodes",
	"conductance_comm",
]


def ego_net_plots(data, graph_set_name, graph_set_params):
	print("Plots for ego-net metrics")
	for ego_metric in ego_metrics:
		make_plot(
			# plot_type=PlotType.bar,
			data=data["ego_net_metrics"].query("metric_name in @ego_metric"),
			graph_filter=graph_set_params['graph_filter'],
			# xlabel="om",
			algo_matches=algo_matches,
			# add_algos=["Ego_PLP_b", "Ego_PLP_e"],
			remove_algo_part=remove_algo_part,
			title="Ego-Net {}".format(ego_metric),
			file_name="ego_partition/metrics/{}_{}".format(ego_metric, graph_set_name),
			one_plot_per_graph=False,
			x=graph_set_params['x'],
			y="value",
			hue="algo",
			plot_args={
				# "style": "metric_name",
			},
			ax_set={
				# "ylim": (0, 1.05),
				"ylim": 0,
			}
		)

	# # Algo parameters on x-axis
	# for ego_metric in ego_metrics:
	# 	break
	# 	x = {"name": "factor", "create_from": "algo", "str_start": "_f-", "str_end": "*"}
	# 	hue = {"name": "exponent", "create_from": "algo", "str_start": "*e-",
	# 	       "str_end": ""}
	# 	make_plot(
	# 		data=data["ego_net_metrics"].query("metric_name in @ego_metric"),
	# 		graphs="LFR_om",
	# 		algo_match=algo_matches,
	# 		title="Ego-Net Metrics, " + ego_metric,
	# 		file_name="ego_partition/2_dim_metrics/" + ego_metric,
	# 		x=x,
	# 		y="value",
	# 		hue=hue,
	# 		plot_args={
	# 		},
	# 		ax_set={
	# 			"ylim": (0, 1.05),
	# 		}
	# 	)


def ego_net_plots_per_graph(data, graph_set_name, graph_set_params):
	# Metrics per Ego-Net
	print("Plots for ego-net metrics per graph")
	for ego_metric in ego_metrics:
		make_plot(
			# plot_type=PlotType.bar,
			data=data["ego_net_ego_metrics"].query("metric_name in @ego_metric"),
			graph_filter=graph_set_params['graph_filter'],
			algo_matches=algo_matches,
			remove_algo_part=remove_algo_part,
			title="Ego-Net {}".format(ego_metric),
			file_name="ego_partition/ego_metrics/{}_{}".format(ego_metric, graph_set_name),
			one_plot_per_graph=True,
			x="ego_net_size",
			y="value",
			hue="algo",
			plot_args={
				# "style": "metric_name",
				# "markersize": 3,
			},
			ax_set={
				# "ylim": (0, 1.05),
				"ylim": 0,
			}
		)
