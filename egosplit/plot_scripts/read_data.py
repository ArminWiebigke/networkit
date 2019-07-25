import glob
from collections import defaultdict

import pandas as pd

from networkit.stopwatch import clockit
from .extract_data_column import create_new_column_from_string


@clockit
def read_data(result_dir):
	if not result_dir[-1] == "/":
		result_dir += "/"
	data = dict()
	for name in ["metrics",
	             "timings",
	             "ego_net_communities",
	             "ego_net_partitions",
	             "ego_net_ego_metrics",
	             "ego_net_metrics",
	             "cover_comm_sizes",
	             "cover_node_comms",
	             "cover_num_comms"
	             ]:
		filter_files = "2019-07-25T20:56*"
		filename = '/{}/{}.result'.format(filter_files, name)
		files = glob.glob(result_dir + filename)
		print(*files)
		# for file in files:
		# 	file_data = pd.read_csv(file, sep="\s+")
		# 	if not name in data:
		# 		data[name] = file_data
		# 	else:
		# 		data[name] = data[name].append(file_data, ignore_index=True)
		if not files:
			print("File " + filename + " not found")
		else:
			data[name] = pd.concat([pd.read_csv(file, sep=",", comment='#') for file in files],
			                       ignore_index=True)

	return data


def create_column_if_missing(data, column):
	if column not in data.columns:
		if column == 'Running Time / Number of Edges':
			data[column] = data['Running Time'] / data['Number of Edges']
		elif column == 'Running Time / Average Degree':
			data[column] = (data['Running Time'] /
			                (data['Number of Edges'] / data['Number of Nodes'])).astype('float64')
		elif column == 'Communities per Node':
			data[column] = (data['om'] - 1) * (data['on'] / data['N']) + 1
		elif column == 'Mixing Factor':
			data[column] = data['mu']
		else:
			create_new_column_from_string(data, column, get_new_column_params()[column])


def get_new_column_params():
	new_columns = {
		# "Communities per Node": {
		# 	"create_from": "Graph Name",
		# 	"str_start": "_om_",
		# 	"str_end": "_",
		# 	'type': float
		# },
		# "Mixing Factor": {
		# 	"create_from": "Graph Name",
		# 	"str_start": "_mu_",
		# 	"str_end": "_",
		# 	'type': float
		# }
	}
	return new_columns
