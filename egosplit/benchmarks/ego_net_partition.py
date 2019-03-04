from networkit.community import EgoSplitting, CoverF1Similarity, PLM, PLP, \
	LPPotts, OLP
from networkit.structures import Partition
from networkit.components import ConnectedComponents
from networkit import none
from egosplit.benchmarks.algorithms import *
from egosplit.benchmarks.execution import *


def bench_ego_net_partitions():
	append = False
	open_mode = 'w'
	if append:
		open_mode = 'a'
	out_file_comm = open('../results/ego_net_partitions.result', open_mode)
	out_file_part = open('../results/ego_net_partition_composition.result', open_mode)
	if not append:
		out_file_comm.write("algo graph comm_size wrong_nodes num_partitions\n")
		out_file_part.write("algo graph partition_size wrong_nodes wrong_percentage\n")

	om = 3
	lfr_params = {
		'N': 2000, 'k': 18 * om, 'maxk': 120, 'minc': 60, 'maxc': 100,
		't1': 2, 't2': 2, 'mu': 0.2, 'on': 2000, 'om': om}
	graph_name = "LFR_om_3"
	graphs = []
	graphs.append(LFRGraph(graph_name, lfr_params))

	partition_algos = OrderedDict()
	partition_algos['PLP'] = [lambda g: PLP(g, 1, 20).run().getPartition()]
	partition_algos['PLM'] = [lambda g: PLM(g, False, 1.0, "none").run().getPartition()]
	partition_algos['LPPotts'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition()]
	partition_algos['Infomap'] = [lambda g: clusterInfomap(g)]

	for algo_name in partition_algos:
		for graph in graphs:
			run_algo(algo_name, partition_algos[algo_name], graph, out_file_comm,
			         out_file_part)


def run_algo(algo_name, partition_algos, graph, out_file_comm, out_file_part):
	graph_name = graph.name
	ground_truth = graph.ground_truth
	ego_algo = EgoSplitting(graph.graph, *partition_algos)
	ego_algo.run()
	ego_net_partitions = ego_algo.getEgoNetPartitions()

	for u in graph.graph.nodes():
		truth_communities = ground_truth.subsetsOf(u)
		partitions, best_communities, community_sizes = calculate_partition_properties(
			ground_truth, truth_communities, ego_net_partitions[u])

		# Are communities split into multiple partitions?
		result_list = check_community_partitioning(
			ground_truth, partitions, best_communities, community_sizes,
			truth_communities)
		for result in result_list:
			out_file_comm.write(
				algo_name + " " + graph_name + " " + str(result["comm_size"])
				+ " " + str(result["wrong_nodes"]) + " " + str(result["num_partitions"])
				+ "\n"
			)

		# Do partitions consist of multiple communities?
		result_list = check_partition_composition(
			ground_truth, partitions, best_communities, community_sizes,
			truth_communities)
		for result in result_list:
			out_file_part.write(
				algo_name + " " + graph_name + " " + str(result["partition_size"])
				+ " " + str(result["wrong_nodes"])
				+ " " + str(result["wrong_percentage"]) + "\n"
			)


# Returns partition->node list, best community for each partition, community->size map
def calculate_partition_properties(ground_truth, truth_communities, partition_map):
	num_partitions = partition_map[none]
	partitions = []
	for i in range(0, num_partitions):
		partitions.append([])
	for v, partition in partition_map.items():
		if v != none:
			partitions[partition].append(v)
	# Remove partitions of size one
	# partitions = [p for p in partitions if len(p) > 1]

	max_community_cnts = {}
	for c in truth_communities:
		max_community_cnts[c] = 0
	community_sizes = {}
	for c in truth_communities:
		community_sizes[c] = 0
	best_communities = [] # for each partition: the ground truth community with the most nodes
	community_cnts = [] # for each partition: number of nodes for each ground truth community
	for i in range(0, len(partitions)):
		partition = partitions[i]
		community_cnts.append(dict())
		for c in truth_communities:
			community_cnts[i][c] = 0
		for v in partition:
			comms = ground_truth.subsetsOf(v)
			for c in comms:
				if c in truth_communities:
					community_cnts[i][c] += 1
		for c in community_cnts[i]:
			max_community_cnts[c] = max(max_community_cnts[c], community_cnts[i][c])
			community_sizes[c] += community_cnts[i][c]

		best_communities.append(dominating_community(community_cnts[i], len(partition)))
	return partitions, best_communities, community_sizes


# For each partition: Get the community with the most nodes. Count the number of nodes
# that are not in that community.
def check_partition_composition(ground_truth, partitions, best_communities,
                                community_sizes, truth_communities):
	result_list = []
	for i in range(0, len(partitions)):
		p = partitions[i]
		best_community = best_communities[i]
		node_cnt = 0
		for u in p:
			if best_community not in ground_truth.subsetsOf(u):
				node_cnt += 1
		result_list.append({
			'partition_size': len(p),
			'wrong_nodes': node_cnt,
			'wrong_percentage': float(node_cnt) / len(p)
		})
		if node_cnt > 0.8 * len(p) and len(p) > 20:
			print("")
			print(p)
			for u in p:
				print([c for c in ground_truth.subsetsOf(u) if c in truth_communities])
			print(best_community)
		return result_list
	return result_list


# Count the number of nodes of the given community in all partitions except the partition
# with the largest number of nodes from that community.
# Don't count nodes that are in multiple communities and that are assigned to one of its
# communities correctly.
def check_community_partitioning(ground_truth, partitions, best_communities,
                                 community_sizes, truth_communities):
	result_list = []
	for c in truth_communities:
		best_partition, best_partition_cnt = find_best_partition(ground_truth, c, partitions)
		node_cnt = 0
		for i in range(0, len(partitions)):
			if i != best_partition:
				for u in partitions[i]:
					comms = ground_truth.subsetsOf(u)
					if c in comms:
						if best_communities[i] not in comms or best_communities[i] == c:
							node_cnt += 1
		result_list.append({
			'comm_size': community_sizes[c],
			'wrong_nodes': node_cnt,
			'num_partitions': best_communities.count(c)
		})
		# if best_communities.count(c) > 2:
		# 	print("Community: ", str(c))
		# 	print("Wrong:     ", str(wrong_nodes[c]) + "/" + str(total_community_cnts[c]))
		# 	print("Partitions:", str(best_communities.count(c)))
		# 	print("Max cnt:   ", max_community_cnts)
		# 	print("Total cnt: ", total_community_cnts)
		# 	print("Best:      ", best_communities)
		# 	print("Counts:    ", community_cnts)
		# 	print("")
	return result_list


# Find the partition with the largest number of nodes that are in the given community
def find_best_partition(ground_truth, community, partitions):
	best_partition = -1
	best_partition_cnt = 0
	for i in range(0, len(partitions)):
		num_nodes = 0
		for u in partitions[i]:
			if community in ground_truth.subsetsOf(u):
				num_nodes += 1
		if num_nodes > best_partition_cnt:
			best_partition_cnt = num_nodes
			best_partition = i
	return best_partition, best_partition_cnt


# For a given partition, return the community with the highest number of nodes.
# If no community has more than 50% of the nodes, return -1.
def dominating_community(community_cnts, num_nodes):
	max_comm = -1
	max_comm_cnt = 0
	for comm, comm_cnt in community_cnts.items():
		if comm_cnt > max_comm_cnt: # and comm_cnt > 0.5 * num_nodes:
			max_comm = comm
			max_comm_cnt = comm_cnt
	return max_comm


bench_ego_net_partitions()