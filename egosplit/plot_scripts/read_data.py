import glob
from collections import defaultdict

import pandas as pd
import os

from networkit.stopwatch import clockit
from egosplit.plot_scripts.extract_data_column import create_new_column_from_string


# @clockit
# def read_data(result_dir, filter_files):
# 	data = defaultdict(lambda f: read_file(filter_files, f, result_dir))
	# for name in ["metrics",
	#              "timings",
	#              "ego_net_communities",
	#              "ego_net_partitions",
	#              "ego_net_ego_metrics",
	#              "ego_net_metrics",
	#              "cover_comm_sizes",
	#              "cover_node_comms",
	#              "cover_num_comms"
	#              ]:
	# 	read_file(filter_files, name, result_dir)

	# return data


class DataReader:
	def __init__(self, result_dir):
		self.data = dict()
		self.result_dir = result_dir

	def __getitem__(self, name):
		if name not in self.data:
			self.data[name] = self.read_file(name)
		return self.data[name]

	def read_file(self, name):
		filename = '{}.result'.format(name)
		file_path = os.path.join(self.result_dir, filename)
		files = glob.glob(file_path)
		if not files:
			raise FileNotFoundError("File " + file_path + " not found")
		else:
			print("Read files:")
			print(*files)
			data = pd.concat([pd.read_csv(file, sep=",", comment='#') for file in files],
			                 ignore_index=True)
		return data


def create_column_if_missing(data, column):
	if column not in data.columns:
		if column == 'Running Time / Number of Edges':
			data[column] = data['Running Time'] / data['Number of Edges']
		elif column == 'Running Time / Average Degree':
			data[column] = (data['Running Time'] /
			                (data['Number of Edges'] / data['Number of Nodes'])).astype('float64')
		elif column == 'time / (n + m)':
			data[column] = 1e6 * (data['Running Time'] /
			                (data['Number of Edges'] + data['Number of Nodes'])).astype('float64')
		elif column == 'Running Time / (n + m)':
			data[column] = 1000 * (data['Running Time'] /
			                (data['Number of Edges'] + data['Number of Nodes'])).astype('float64')
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
