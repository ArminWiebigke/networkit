from collections import OrderedDict
from egosplit.external import getFacebookGraph, getAmazonGraph, getDBLPGraph
from egosplit.benchmarks.data_structures.graph import LFRGraph, ReadGraph, BenchGraph


def get_graphs(graph_sets):
	return GraphSetsConfig.get_graphs(graph_sets)


class GraphSetsConfig:
	@staticmethod
	def get_graphs(graph_sets):
		all_graphs = []
		for graph_set in graph_sets:
			graphs = []
			if graph_set == "om":
				graphs = create_om_graphs(1, 7)
			elif graph_set == "mu":
				graphs = create_mu_graphs(10, 70, 10)
			elif graph_set == "overlap":
				graphs = create_overlap_graphs(20, 80)
			elif graph_set == "test":
				graphs = create_om_graphs(3, 3)
			elif graph_set == 'facebook':
				graphs = facebook_graphs()

			all_graphs.extend(graphs)
		return all_graphs


def om_graphs():
	return create_om_graphs(2, 5)


def mu_graphs():
	return create_mu_graphs(20, 40, 5)


def overlap_graphs():
	return create_overlap_graphs()


def test_graphs():
	return create_om_graphs(3, 3)


def facebook_graphs():
	graphs = []
	graphs.append(
		ReadGraph(lambda: getFacebookGraph('Caltech36', clean=True), 'FB_1_Caltech36'))  # 769 nodes
	graphs.append(
		ReadGraph(lambda: getFacebookGraph('Smith60', clean=True), 'FB_2_Smith60'))  # 3k nodes
	graphs.append(
		ReadGraph(lambda: getFacebookGraph('Rice31', clean=True), 'FB_3_Rice31'))  # 4k nodes
	graphs.append(
		ReadGraph(lambda: getFacebookGraph('Auburn71', clean=True), 'FB_4_Auburn71'))  # 18k nodes
	# graphs.append(ReadGraph(lambda: getFacebookGraph('Oklahoma97', clean=True), 'FB_Oklahoma97'))  # 17k nodes
	return graphs


def large_graphs():
	graphs = []
	graphs.append(ReadGraph(lambda: getAmazonGraph(), 'Amazon'))
	graphs.append(ReadGraph(lambda: getDBLPGraph(), 'DBLP'))
	# graphs.append(ReadGraph(lambda: getLiveJournalGraph(), 'LiveJournal'))
	# graphs.append(ReadGraph(lambda: getOrkutGraph(), 'Orkut'))
	return graphs


def create_mu_graphs(min_mu=10, max_mu=80, mu_step=5):
	""" Scale the mixing factor. """
	lfr_graph_args = OrderedDict()
	for mu_factor in range(min_mu, max_mu + 1, mu_step):
		om = 3
		mu = 0.01 * mu_factor
		k = 10 * om  # Number of neighbors per community independent of Mixing Factor
		k /= (1 - mu)
		maxk = 20 + 10 * om  # Scale max degree with average degree
		maxk /= (1 - mu)
		name = 'mu_{:02.0f}'.format(mu_factor)
		lfr_graph_args[name] = {
			'N': 2000, 'k': k, 'maxk': maxk, 'minc': 30, 'maxc': 60,
			't1': 2, 't2': 2, 'mu': 0.01 * mu_factor, 'on': 2000, 'om': om}

	graphs = create_lfr_graphs(lfr_graph_args)
	return graphs


def create_om_graphs(min_om=1, max_om=7):
	""" Scale the number of communities per node. """
	lfr_graph_args = OrderedDict()
	N = 2000
	for om in range(min_om, max_om + 1):
		on = N
		name = 'om_{}'.format(om)
		graph_name, graph_args = get_om_parameter(N, on, om, name)
		lfr_graph_args[graph_name] = graph_args

	graphs = create_lfr_graphs(lfr_graph_args)
	return graphs


def get_om_parameter(N, on, om, name):
	mu = 0.25
	minc = 30
	maxc = 60
	pcnt_overlap = on / N
	avg_comms = pcnt_overlap * om + (1 - pcnt_overlap)
	k = 10 * avg_comms  # Number of neighbors per community independent of Mixing Factor
	k /= (1 - mu)
	maxk = 20 + 10 * avg_comms  # 2 * k  # Scale max degree with average degree
	maxk /= (1 - mu)
	return (name, {
		'N': N, 'k': k, 'maxk': maxk, 'minc': minc, 'maxc': maxc,
		't1': 2, 't2': 2, 'mu': mu, 'on': on, 'om': om})


def create_overlap_graphs(min_overlap=20, max_overlap=80, step=20):
	lfr_graph_args = OrderedDict()
	N = 2000
	for overlap in range(min_overlap, max_overlap + 1, step):
		on = N * overlap / 100
		avg_comms = 1 + overlap / 100
		om = 2
		name = 'om_{:.2f}'.format(avg_comms)
		graph_name, graph_args = get_om_parameter(N, on, om, name)
		lfr_graph_args[graph_name] = graph_args

	graphs = create_lfr_graphs(lfr_graph_args)
	return graphs


def create_lfr_graphs(graphArgs):
	graphs = []
	for argsName in graphArgs.keys():
		args = graphArgs[argsName]
		graph_name = 'LFR_' + argsName
		print(args)
		graph_wrapper = LFRGraph(graph_name, args)
		graphs.append(graph_wrapper)
	return graphs
