import numpy as np

from .config import *
from .utils import filter_data


def make_plot(data, graphs, algo_match=None, xlabel=None, remove_algo_part="", file_name="", title="",
              plot_args=None, ax_set=None, plot_type="line", add_algos=None, x=None, hue=None, y=None):
	filtered_data = data.query("graph.str.contains(@graphs)").copy()
	if algo_match is not None:
		algo_list = list(filtered_data.query("algo.str.contains(@algo_match)").groupby("algo").mean()
	                     .index.values)
	else:
		algo_list = []
	if add_algos:
		algo_list.extend(add_algos)
	filtered_data.query("algo in @algo_list", inplace=True)

	if len(filtered_data) is 0:
		return

	new_columns = []
	if type(x) != str:
		new_columns.append(x)
	if type(hue) != str:
		new_columns.append(hue)
	for new_column in new_columns:
		for row in filtered_data.itertuples(index=True, name='Pandas'):
			index = row.Index
			from_string = filtered_data.loc[index, new_column["create_from"]]
			start_idx = from_string.find(new_column["str_start"]) + len(new_column["str_start"])
			if new_column["str_end"] == "":
				end_idx = len(from_string)
			else:
				end_idx = from_string.find(new_column["str_end"], start_idx)
			value = from_string[start_idx:end_idx]
			try:
				value = float(value)
			except ValueError:
				pass
			filtered_data.loc[index, new_column["name"]] = value

	if type(x) != str:
		x = x["name"]
	if type(hue) != str:
		hue = hue["name"]

	if plot_type == "line":
		default_plot_args = {
			"markers": True,
			"linewidth": 2,
			"ci": None,
			"style": hue,
		}
	elif plot_type == "swarm":
		default_plot_args = {
			"size": 2,
		}
	else:
		raise RuntimeError("No valid plot type!")
	plot_args = plot_args or {}
	plot_args = {
		**default_plot_args,
		"x": x,
		"y": y,
		"hue": hue,
		**plot_args,
	}

	# Make one plot per graph if graph is not on the x-axis
	if plot_args["x"] == "graph":
		graph_list = [graphs]
	else:
		graph_list = filtered_data.groupby("graph").mean().index.values

	for graph in graph_list:
		if len(graph_list) == 1:
			graph_data = filtered_data
		else:
			graph_data = filtered_data.query("graph == @graph")

		fig, ax = plt.subplots()
		if plot_type == "line":
			sns.lineplot(
				**plot_args,
				data=graph_data,
				ax=ax,
			)
		elif plot_type == "swarm":
			sns.swarmplot(
				**plot_args,
				data=graph_data,
				ax=ax,
			)
		else:
			raise RuntimeError("No valid plot type!")

		sns.despine(ax=ax)
		ax_set = ax_set or {}
		xlabel = xlabel or ax.get_xlabel()
		ax_set = {
			"xlabel": xlabel,
			**ax_set,
		}
		ax.set(
			**ax_set,
		)
		set_xticklabels(ax, xlabel)
		fig.suptitle(title + ", " + graph)
		legend_handles, legend_labels = ax.get_legend_handles_labels()
		if type(algo_match) == str and remove_algo_part == "":
			remove_algo_part = algo_match
		legend_labels = [l.replace(remove_algo_part, '') for l in legend_labels]
		set_layout(ax, legend_handles, legend_labels)
		fig.savefig(file_prefix + file_name + "_" + graph + ".pdf")
		plt.close(fig)


def two_dim_algo_plot(data, graphs, algos, remove_algo_part="", file_name="", title="",
                      plot_args=None, ax_set=None, plot_type="line", add_algos=None,
                      x=None, hue=None):
	filtered_data = filter_data(data, graphs, algos)
	if len(filtered_data) is 0:
		return

	add_algos = add_algos or []
	for algo in add_algos:
		filtered_data = filtered_data.append(filter_data(data, graphs, [algo]))

	for new_column in [x, hue]:
		for row in filtered_data.itertuples(index=True, name='Pandas'):
			index = row.Index
			from_string = filtered_data.loc[index, new_column["create_from"]]
			start_idx = from_string.find(new_column["str_start"]) + len(new_column["str_start"])
			if new_column["str_end"] == "":
				end_idx = len(from_string)
			else:
				end_idx = from_string.find(new_column["str_end"], start_idx)
			value = from_string[start_idx:end_idx]
			filtered_data.loc[index, new_column["name"]] = float(value)

	graph_list = filtered_data.groupby("graph").mean().index.values

	default_plot_args = {
		"markers": True,
		"linewidth": 2,
		"ci": None,
	}

	plot_args = plot_args or {}
	plot_args = {
		**default_plot_args,
		"x": x["name"],
		"hue": hue["name"],
		"style": hue["name"],
		**plot_args,
	}

	for graph in graph_list:
		graph_data = filtered_data.query("graph == @graph")
		num_hues = len(filtered_data.groupby(hue["name"]).mean().index.values)
		fig, ax = plt.subplots()
		sns.lineplot(
			**plot_args,
			data=graph_data,
			ax=ax,
			palette=sns.light_palette("green", num_hues),
		)

		sns.despine(ax=ax)
		ax_set = ax_set or {}
		ax.set(
			**ax_set,
		)
		fig.suptitle(title + ", " + graph)
		set_layout(ax)
		fig.savefig(file_prefix + file_name + "_" + graph + ".pdf")
		plt.close(fig)
