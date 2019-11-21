from egosplit.benchmarks.execution.algo_creation import EgoSplitClusteringAlgorithmsConfig, \
	EgoSplitParameterConfig, OtherAlgorithms
from egosplit.benchmarks.execution.cleanup import CleanUpConfig
from egosplit.benchmarks.execution.graph_creation import GraphSetsConfig
from egosplit.benchmarks.plot_scripts.bench_config import PlotGraphSetConfig, PlotAlgoSetConfig
from egosplit.benchmarks.plot_scripts.create_plots import PlotSetConfig


class BenchmarkSet:
	config = None

	def __init__(self):
		self.name = None  # Name of the configuration
		# Benchmark configuration
		self.result_dir = None  # suffix of the direction where the results are stored
		self.plot_dir = None  # the subdirectory in the plots folder where plots will be created
		self.store_ego_nets = False  # if True, store the ego-nets and calculate metrics for them
		self.score_per_egonet = False  # if True, store the metrics for each ego-net, else only the average over all ego-nets is stored
		self[EgoSplitClusteringAlgorithmsConfig] = None  # clustering algorithms for EgoSplitting
		self[EgoSplitParameterConfig] = None  # parameters of EgoSplitting
		self[GraphSetsConfig] = None  # input graphs used in the benchmarks
		self[CleanUpConfig] = 'No Cleanup'  # clean up the result of EgoSplitting
		self.stream_to_gephi = False  # stream the results to Gephi for visualization
		self[OtherAlgorithms] = None  # additionally executed algorithms (not EgoSplitting)
		self.calc_f1_per_comm = False
		self.iterations = 10
		self.time_limit = 3600  # time limit in seconds
		# Plot configuration
		self.remove_algo_parts = []  # remove parts of the algorithm name when displayed in the plots
		self.replace_legend = {}  # Replace legend entries in the plots
		self.hue = 'Algorithm'  # the hue used in the plots
		self[PlotGraphSetConfig] = None  # graph sets used in the plots,
		self[PlotSetConfig] = None  # which plots to create
		self[PlotAlgoSetConfig] = None  # algorithm sets for the plots

		# Set config from derived class, prevent any new attributes that are not defined above
		self._frozen = True
		for key, value in self.config.items():
			self[key] = value

	def __setitem__(self, key, value):
		self.__setattr__(key, value)

	def __getitem__(self, item):
		return self.__dict__[item]

	def __setattr__(self, key, value):
		if hasattr(self, '_frozen') and not key in self.__dict__.keys():
			raise AttributeError('BenchmarkSet has no attribute \'{}\' (set by {})'.format(
				key, self.__class__.__name__))
		self.__dict__[key] = value
