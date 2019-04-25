def filter_data(data, graphs, algos):
	filtered_data = data.query("graph.str.contains(@graphs)").copy()
	if type(algos) == list:
		filtered_data.query("algo in @algos", inplace=True)
	else:
		filtered_data.query("algo.str.contains(@algos)", inplace=True)
	return filtered_data
