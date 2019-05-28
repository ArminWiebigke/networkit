import matplotlib.pyplot as plt
import seaborn as sns

from enum import Enum

from .config import set_xticklabels, set_layout, file_prefix


class PlotType(Enum):
	line = 1
	bar = 2
	swarm = 3
	violin = 4


# Create a plot
def make_plot(data, graphs, algo_match="", xlabel=None, remove_algo_part="",
              file_name="", title="",
              plot_args=None, ax_set=None, plot_type=PlotType.line, add_algos=None, x=None,
              hue=None, y=None, one_plot_per_graph=False):
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
	gt = "Ground_Truth"
	if gt in algo_set:
		algo_set.remove(gt)
		algo_list = sorted(list(algo_set)) + [gt]
	else:
		algo_list = sorted(list(algo_set))
	filtered_data.query("algo in @algo_list", inplace=True)
	if len(filtered_data) is 0:
		return

	# Create new data columns for x or hue if necessary
	new_columns = []
	if type(x) != str and x:
		new_columns.append(x)
	if type(hue) != str and hue:
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
	if type(x) != str and x:
		x = x["name"]
	if type(hue) != str and hue:
		hue = hue["name"]

	# Set plot args
	if plot_type == PlotType.line:
		default_plot_args = {
			"markers": True,
			"linewidth": 2,
			"ci": None,
			"style": hue,
		}
	elif plot_type == PlotType.swarm:
		default_plot_args = {
			"size": 2,
		}
	elif plot_type is PlotType.bar:
		default_plot_args = {
		}
	elif plot_type is PlotType.violin:
		default_plot_args = {
			"scale": "count",
			"linewidth": 1,
			"inner": None,
			"cut": 0,
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
	if plot_args["hue"] == "algo":
		plot_args = {"hue_order": algo_list, **plot_args}
		if plot_type is PlotType.line:
			plot_args = {"style_order": algo_list, **plot_args}
		elif plot_args["x"] == "algo" and plot_type is PlotType.swarm:
			plot_args = {"order": algo_list, **plot_args}

	# Make one plot per graph if graph is not on the x-axis
	if plot_args["x"] != "graph" or one_plot_per_graph:
		graph_list = filtered_data.groupby("graph").mean().index.values
	else:
		graph_list = [graphs]

	# Create plots
	for graph in graph_list:
		if len(graph_list) == 1:
			graph_data = filtered_data
		else:
			graph_data = filtered_data.query("graph == @graph")

		fig, ax = plt.subplots()
		plot_args = {
			**plot_args,
			"data": graph_data,
			"ax": ax,
		}
		# Draw correct plot type
		plot_functions = {
			PlotType.line: sns.lineplot,
			PlotType.bar: sns.barplot,
			PlotType.swarm: sns.swarmplot,
			PlotType.violin: sns.violinplot,
		}
		plot_functions[plot_type](**plot_args)

		sns.despine(ax=ax)
		ax_set = ax_set or {}
		old_xlabel = ax.get_xlabel()
		xlabel = xlabel or old_xlabel
		ax_set = {
			"xlabel": xlabel,
			**ax_set,
		}
		ax.set(
			**ax_set,
		)
		fig.canvas.draw()
		if xlabel != old_xlabel:
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
