import sys

from networkit import gephi

from ..algorithms import EgoSplitAlgorithm


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


def mark_direct_neighbors(graph, egonet, node_id):
	neighbors = graph.neighbors(node_id)
	neighbor_map = dict()
	for u in egonet.nodes():
		if u in neighbors:
			neighbor_map[u] = 1
		else:
			neighbor_map[u] = 0
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

				neighbors = mark_direct_neighbors(graph.graph, egonet, node_id)
				client.exportNodeValues(egonet, neighbors, "neighbors")

				partition = benchmark.algo.ego_net_partition_of(node_id)
				client.exportNodeValues(egonet, partition, benchmark.algo.name)

			command = input("Press Enter to proceed, or q + Enter to quit")
			if command == "q":
				keep_streaming = False
				break