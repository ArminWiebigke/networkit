from plot_scripts.comm_sizes import comm_sizes_ego,	comm_sizes_lfr
from plot_scripts.comms_per_node import node_comms_ego, node_comms_lfr
from plot_scripts.egonet import egonet_comm_partition, egonet_partition_composition
from plot_scripts.metrics import metrics_lfr
from plot_scripts.num_comms import num_comms_lfr
from plot_scripts.partition_counts import partition_counts
from plot_scripts.read_data import read_data
from plot_scripts.config import set_sns_style

data, metrics, all_graph_names, all_algo_names = read_data()
set_sns_style()

# Create plots
metrics_lfr(data)

# num_comms_lfr(data)
#
# comm_sizes_lfr(data)
# comm_sizes_ego(data)
#
# node_comms_lfr(data)
# node_comms_ego(data)
#
# partition_counts(data)
#
# egonet_comm_partition(data, "LFR_om_3")
# egonet_partition_composition(data, "LFR_om_3")

# plt.show()


# def filter_data(data, conditions):
# 	for cond in conditions:
# 		data = data[data[cond[0]] == cond[1]]
# 	return data
#
#
# def get_value(base_data, x_name, x_val, y_name):
# 	graph_data = filter_data(base_data, [(x_name, x_val)])
# 	grouped_data = graph_data.groupby(x_name)
# 	y = grouped_data.mean()[y_name].values
# 	val = 0.0
# 	if len(y) == 1:
# 		val = y[0]
# 	elif len(y) != 0:
# 		print("ERROR", y, len(y))
# 	return val
#
#
# def metrics_real():
# 	file_append = 'real'
# 	filtered_data = data["metrics"].query("graph.str.contains('real')")
# 	for metric in metrics:
# 		fig, ax = plt.subplots()
# 		sns.barplot(x="graph",
# 					y=metric,
# 					hue="algo",
# 					palette="bright",
# 					data=filtered_data,
# 					ax=ax)
# 		sns.despine(left=True, ax=ax)
# 		ax.set(
# 			ylabel=metric_names[metric]["y_val"],
# 			xticklabels=[label.get_text()[5:] for label in ax.get_xticklabels()]
# 		)
# 		if metric is "time":
# 			ax.set_yscale('log')
# 		ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3)
# 		fig.suptitle(metric_names[metric]["description"] + " real graphs")
# 		plt.tight_layout(rect=(0, 0, 1, 0.96))
# 		fig.savefig(file_prefix + "metrics/" + metric_names[metric]["file_name"] + "_" + file_append + ".pdf")
