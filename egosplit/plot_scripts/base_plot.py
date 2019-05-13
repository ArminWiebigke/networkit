import numpy as np

from .config import *
from .utils import filter_data


def make_plot(data, graphs, algo_match=None, xlabel=None, remove_algo_part="",
              file_name="", title="",
              plot_args=None, ax_set=None, plot_type="line", add_algos=None, x=None,
              hue=None, y=None):
	# Filter data
	filtered_data = data.query("graph.str.contains(@graphs)").copy()
	if algo_match is not None:
		f_data = filtered_data.query("algo.str.contains(@algo_match)")
		algo_set = set(f_data.groupby("algo").mean().index.values)
	else:
		algo_set = set()
	if add_algos:
		for algo in add_algos:
			algo_set.add(algo)
	algo_list = sorted(list(algo_set))
	filtered_data.query("algo in @algo_list", inplace=True)
	if len(filtered_data) is 0:
		return

	# Create new columns for x or hue if necessary
	new_columns = []
	if type(x) != str:
		new_columns.append(x)
	if type(hue) != str:
		new_columns.append(hue)
	for new_column in new_columns:
		for row in filtered_data.itertuples(index=True, name='Pandas'):
			index = row.Index
			from_string = filtered_data.loc[index, new_column["create_from"]]
			start_idx = from_string.find(new_column["str_start"]) + len(
				new_column["str_start"])
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

	# Set plot args
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
		# "sort": True,
		**plot_args,
	}
	if plot_args["hue"] == "algo":
		plot_args = {"hue_order": algo_list, **plot_args}
		if plot_type == "line":
			plot_args = {"style_order": algo_list, **plot_args}
		if plot_type == "swarm":
			plot_args = {"order": algo_list, **plot_args}

	# Make one plot per graph if graph is not on the x-axis
	if plot_args["x"] == "graph":
		graph_list = [graphs]
	else:
		graph_list = filtered_data.groupby("graph").mean().index.values

	# Create plots
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
		fig.canvas.draw()
		set_xticklabels(ax, xlabel + "_")
		fig.suptitle(title + ", " + graph)
		legend_handles, legend_labels = ax.get_legend_handles_labels()
		if type(algo_match) == str and remove_algo_part == "":
			remove_algo_part = algo_match
		if isinstance(remove_algo_part, str):
			remove_algo_part = [remove_algo_part]
		for remove in remove_algo_part:
			legend_labels = [l.replace(remove, '') for l in legend_labels]
		set_layout(ax, legend_handles, legend_labels)
		fig.savefig(file_prefix + file_name + "_" + graph + ".pdf")
		plt.close(fig)
