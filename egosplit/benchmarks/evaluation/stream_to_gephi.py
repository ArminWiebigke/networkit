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
	for graph in graphs:
		G = graph.graph
		for _ in range(1):
			node_id = G.randomNode()

			for benchmark in benchmarks:
				if benchmark.graph is not graph \
						or not isinstance(benchmark.algo, EgoSplitAlgorithm):
					continue
				egonet = benchmark.algo.getEgoNet(node_id)
				if len(egonet.nodes()) == 0:
					print("No egonet stored")
					continue
				egonet.indexEdges()

				client = gephi.streaming.GephiStreamingClient(
					url='http://localhost:8080/workspace' + str(workspace_id))
				client.exportGraph(egonet)

				assert egonet.hasEdgeIds()
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

				partition = benchmark.algo.getEgoNetPartitions()[node_id]
				client.exportNodeValues(egonet, partition, benchmark.algo.name)

				workspace_id += 1
