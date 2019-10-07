from collections import defaultdict

from networkit import gephi

from egosplit.benchmarks.data_structures.algorithms import EgoSplitAlgorithm


def create_ground_truth_partition(node_id, egonet, ground_truth):
	gt_comms = ground_truth.subsetsOf(node_id)
	partition = dict()
	for u in egonet.nodes():
		comms = ground_truth.subsetsOf(u).intersection(gt_comms)
		if len(comms) == 0:
			color = 1000
		elif len(comms) > 1:
			color = 1001
		else:
			color = list(comms)[0]
		partition[u] = color
	return partition


def create_ground_truth_coloring(node_id, egonet, ground_truth, color_dict):
	gt_comms = ground_truth.subsetsOf(node_id)
	partition = dict()
	for u in egonet.nodes():
		comms = ground_truth.subsetsOf(u).intersection(gt_comms)
		if len(comms) == 0:
			color = color_dict['none']
		elif len(comms) > 1:
			color = color_dict['multiple']
		else:
			comm = comms.pop()
			# if comm not in color_dict:
			# 	color_dict[comm] = colors.pop(0)
			color = color_dict[comm]

		partition[u] = color
	return partition


def mark_direct_neighbors(graph, egonet, node_id):
	neighbors = graph.neighbors(node_id)
	neighbor_map = dict()
	for u in egonet.nodes():
		if u in neighbors:
			neighbor_map[u] = 20
		else:
			neighbor_map[u] = 10
	return neighbor_map


def extend_graph_to_node(graph, u):
	for _ in range(graph.upperNodeIdBound(), u + 1):
		v = graph.addNode()
		graph.removeNode(v)


def stream_partition(graphs, benchmarks):
	workspace_id = 1
	benchmarks_and_clients = []
	for benchmark in benchmarks:
		if not isinstance(benchmark.algo, EgoSplitAlgorithm):
			continue
		egonet = benchmark.algo.ego_net_of(0)
		if len(egonet.nodes()) == 0:
			print("No egonet stored")
			continue
		print(workspace_id)
		client = gephi.streaming.GephiStreamingClient(
			url='http://localhost:8080/workspace' + str(workspace_id))
		benchmarks_and_clients.append((benchmark, client))
		workspace_id += 1

	keep_streaming = True
	while keep_streaming:
		colors = [
			'#3379D9',  # blue
			'#2AA925',  # green
			'#FF2000',  # red
			'#000000',
			'#000000',
			'#000000',
			'#000000',
		]
		color_dict = defaultdict(lambda: colors.pop(0))
		color_dict['none'] = '#cccccc'
		color_dict['multiple'] = '#333333'

		for graph in graphs:
			node_id = graph.graph.randomNode()
			for benchmark, client in benchmarks_and_clients:
				egonet = benchmark.algo.ego_net_of(node_id)
				egonet.indexEdges()
				client.exportGraph(egonet)

				edge_weights = [0 for _ in range(egonet.upperEdgeIdBound())]

				def set_weight(edge_id, weight):
					edge_weights[edge_id] = weight

				egonet.forEdges(lambda u, v, weight, edge_id: set_weight(edge_id, weight))
				client.exportEdgeValues(egonet, edge_weights, "weight")

				gt_partition = create_ground_truth_partition(node_id, egonet,
				                                             graph.ground_truth)
				client.exportNodeValues(egonet, gt_partition, "ground_truth")

				gt_colors = create_ground_truth_coloring(node_id, egonet, graph.ground_truth,
				                                         color_dict)
				client.exportNodeValues(egonet, gt_colors, "color")
				client.exportNodeValues(egonet, gt_colors, "color2")

				neighbors = mark_direct_neighbors(graph.graph, egonet, node_id)
				client.exportNodeValues(egonet, neighbors, "neighbors")
				client.exportNodeValues(egonet, neighbors, "size")

				partition = benchmark.algo.ego_net_partition_of(node_id)
				client.exportNodeValues(egonet, partition, benchmark.get_algo_name())

			command = input("Press Enter to proceed, or q + Enter to quit")
			if command == "q":
				keep_streaming = False
				break