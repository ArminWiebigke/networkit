import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def read_data():
	result_dir = "results/"
	data = dict()

	for name in ["metrics",
	             "ego_net_communities", "ego_net_partitions", "ego_net_ego_metrics",
	             "ego_net_metrics",
	             "cover_comm_sizes", "cover_node_comms", "cover_num_comms"]:
		filename = name + '.result'
		try:
			data[name] = pd.read_csv(result_dir + filename, sep="\s+")
		except FileNotFoundError:
			print("File " + filename + " not found")

	return data