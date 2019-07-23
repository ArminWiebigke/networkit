import pandas as pd
from .extract_data_column import create_new_column


def read_data(result_dir):
	if not result_dir[-1] == "/":
		result_dir += "/"
	data = dict()
	for name in [
		"metrics",
		"ego_net_communities",
		"ego_net_partitions",
		"ego_net_ego_metrics",
		"ego_net_metrics",
		"cover_comm_sizes",
		"cover_node_comms",
		"cover_num_comms"
	]:
		filename = name + '.result'
		try:
			data[name] = pd.read_csv(result_dir + filename, sep="\s+")
		except FileNotFoundError:
			print("File " + filename + " not found")

	return data


def create_column_if_missing(data, column):
	if not column in data.columns:
		create_new_column(data, column, get_new_column_params()[column])


def get_new_column_params():
	new_columns = {
		"communities_per_node": {
			"create_from": "graph",
			"str_start": "_om_",
			"str_end": "_",
			'type': float
		},
		"mixing_factor": {
			"create_from": "graph",
			"str_start": "_mu_",
			"str_end": "_",
			'type': float
		}
	}
	return new_columns
