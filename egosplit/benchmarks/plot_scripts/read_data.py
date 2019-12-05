import glob

import pandas as pd
import os


class DataReader:
	def __init__(self, result_dir):
		self.data = dict()
		self.result_dir = result_dir

	def __getitem__(self, name):
		if name not in self.data:
			self.data[name] = self.read_file(name)
		return self.data[name]

	def read_file(self, name):
		files = []
		for filename in ['*/{}.result', '{}.result']:
			filename = filename.format(name)
			file_path = os.path.join(self.result_dir, filename)
			files += glob.glob(file_path)

		if not files:
			raise FileNotFoundError('File ' + file_path + ' not found')
		else:
			print('Read files:')
			print(*files)
			data = pd.concat([pd.read_csv(file, sep=',', comment='#') for file in files],
			                 ignore_index=True)
		return data


def create_column_if_missing(data, column):
	if column not in data.columns:
		if column == 'Running Time / Number of Edges':
			data[column] = data['Running Time'] / data['Number of Edges']
		elif column == 'Running Time / Average Degree':
			data[column] = (data['Running Time'] /
			                (data['Number of Edges'] / data['Number of Nodes'])).astype('float64')
		elif column == 'Running Time / (n + m)':
			data[column] = 1e6 * (data['Running Time'] /
			                       (data['Number of Edges'] + data['Number of Nodes'])).astype('float64')
		elif column == 'Communities per Node':
			data[column] = (data['om'] - 1) * (data['on'] / data['N']) + 1
		elif column == 'Mixing Factor':
			data[column] = data['mu']
		else:
			raise ValueError('Unknown column \'{}\''.format(column))
