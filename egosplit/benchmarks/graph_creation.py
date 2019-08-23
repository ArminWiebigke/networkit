from collections import OrderedDict

from egosplit.external import getFacebookGraph, getAmazonGraph, getDBLPGraph, getAmazonGraphAll, \
	getOrkutGraph, getLiveJournalGraph
from egosplit.benchmarks.graph import LFRGraph, ReadGraph


def get_graphs(graph_sets, iterations):
	graphs = []
	for graph_set in graph_sets:
		if graph_set == "om":
			graphs.extend(om_graphs(iterations, 1, 7))
		if graph_set == "test":
			graphs.extend(om_graphs(1, 3, 3))
		if graph_set == "overlap":
			graphs.extend(overlap_graphs(iterations))
		if graph_set == "mu":
			graphs.extend(mu_graphs(iterations))
		if graph_set == "large":
			graphs.extend(large_graphs())
		if graph_set == "facebook":
			graphs.extend(facebook_graphs())
	# for graph in graphs:
	# 	graph.graph.indexEdges()
	return graphs


def facebook_graphs():
	graphs = []
	graphs.append(ReadGraph(lambda: getFacebookGraph('Caltech36', clean=True), 'FB_1_Caltech36'))  # 769 nodes
	graphs.append(ReadGraph(lambda: getFacebookGraph('Smith60', clean=True), 'FB_2_Smith60'))  # 3k nodes
	graphs.append(ReadGraph(lambda: getFacebookGraph('Rice31', clean=True), 'FB_3_Rice31'))  # 4k nodes
	graphs.append(ReadGraph(lambda: getFacebookGraph('Auburn71', clean=True), 'FB_4_Auburn71'))  # 18k nodes
	# graphs.append(ReadGraph(lambda: getFacebookGraph('Oklahoma97', clean=True), 'FB_Oklahoma97'))  # 17k nodes
	return graphs


def large_graphs():
	graphs = []
	graphs.append(ReadGraph(lambda: getAmazonGraph(), 'Amazon'))
	graphs.append(ReadGraph(lambda: getDBLPGraph(), 'DBLP'))
	# graphs.append(ReadGraph(lambda: getLiveJournalGraph(), 'LiveJournal'))
	# graphs.append(ReadGraph(lambda: getOrkutGraph(), 'Orkut'))
	return graphs


def mu_graphs(iterations):
	""" Scale the mixing factor. """
	lfr_graph_args = OrderedDict()
	for mu_factor in range(10, 81, 5):
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

	graphs = create_LFR_graphs(lfr_graph_args, iterations)
	return graphs


def om_graphs(iterations, min_om=1, max_om=7):
	""" Scale the number of communities per node. """
	lfr_graph_args = OrderedDict()
	N = 2000
	for om in range(min_om, max_om + 1):
		on = N
		name = 'om_{}'.format(om)
		graph_name, graph_args = om_graph(N, on, om, name)
		lfr_graph_args[graph_name] = graph_args

	graphs = create_LFR_graphs(lfr_graph_args, iterations)
	return graphs


def overlap_graphs(iterations):
	lfr_graph_args = OrderedDict()
	N = 2000
	for overlap in range(20, 100, 20):
		on = N * overlap / 100
		avg_comms = 1 + overlap / 100
		om = 2
		name = 'om_{:.2f}'.format(avg_comms)
		graph_name, graph_args = om_graph(N, on, om, name)
		lfr_graph_args[graph_name] = graph_args

	graphs = create_LFR_graphs(lfr_graph_args, iterations)
	return graphs


def om_graph(N, on, om, name):
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


def create_LFR_graphs(graphArgs, iterations):
	graphs = []
	for argsName in graphArgs.keys():
		args = graphArgs[argsName]
		graph_name = 'LFR_' + argsName
		print(args)
		for i in range(iterations):
			graph_wrapper = LFRGraph(graph_name, args)
			graphs.append(graph_wrapper)
	return graphs
