from networkit.community import EgoSplitting, CoverF1Similarity, PLM, PLP, \
	LPPotts, OLP
from networkit.structures import Partition
from networkit.components import ConnectedComponents
from networkit import none
from egosplit.benchmarks.algorithms import *
from egosplit.benchmarks.output import *


def analyse_ego_net_partitions(benchmarks, result_dir, append):
	open_mode = 'w'
	if append:
		open_mode = 'a'
	out_file_comm = open(result_dir + 'ego_net_partitions.result', open_mode)
	out_file_part = open(result_dir + 'ego_net_partition_composition.result', open_mode)
	if not append:
		print_headers(out_file_comm, out_file_part)

	for benchmark in benchmarks:
		analyse_ego_net_partition(benchmark, out_file_comm, out_file_part)


def print_headers(out_file_comm, out_file_part):
	out_file_comm.write(create_line(
		"algo",
		"graph",
		"comm_size",
		"wrong_nodes",
		"wrong_percentage",
		"num_partitions"
	))
	out_file_part.write(create_line(
		"algo",
		"graph",
		"partition_size",
		"wrong_nodes",
		"wrong_nodes_gt",
		"wrong_nodes_other",
		"wrong_percentage",
		"wrong_percentage_gt",
		"wrong_percentage_other"
	))


def analyse_ego_net_partition(benchmark, out_file_comm, out_file_part):
	try:
		ego_net_partitions = benchmark.algo.getEgoNetPartitions()
	except AttributeError:
		return
	ground_truth = benchmark.graph.ground_truth
	algo_name = benchmark.algo.name
	graph_name = benchmark.graph.name
	for u in benchmark.graph.graph.nodes():
		truth_communities = set(ground_truth.subsetsOf(u))
		partitions, best_communities, community_sizes = calculate_partition_properties(
			ground_truth, truth_communities, ego_net_partitions[u])

		# Are communities split into multiple partitions?
		result_list = check_community_partitioning(
			ground_truth, partitions, best_communities, community_sizes,
			truth_communities)
		for result in result_list:
			out_file_comm.write(create_line(
				algo_name,
				graph_name,
				result["comm_size"],
				result["wrong_nodes"],
				result["wrong_percentage"],
				result["num_partitions"]
			))

		# Do partitions consist of multiple communities?
		result_list = check_partition_composition(
			ground_truth, partitions, best_communities, truth_communities)
		for result in result_list:
			out_file_part.write(create_line(
				algo_name,
				graph_name,
				result["partition_size"],
				result["wrong_nodes"],
				result["wrong_nodes_gt"],
				result["wrong_nodes_other"],
				result["wrong_percentage"],
				result["wrong_percentage_gt"],
				result["wrong_percentage_other"]
			))


# Returns partition->node list, best community for each partition, community->size map
def calculate_partition_properties(ground_truth, truth_communities, partition_map):
	max_partition = 0
	for id in partition_map.values():
		max_partition = max(max_partition, id)
	partitions = [[] for _ in range(max_partition + 1)]
	for v, p_id in partition_map.items():
		if v != none and p_id != none: # TODO: Why does leidenalg assign partition = none, but no other algorithm?
			partitions[p_id].append(v)
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
	for i in range(len(partitions)):
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
                                truth_communities):
	result_list = []
	for i in range(len(partitions)):
		p = partitions[i]
		num_partitions = len(p)
		best_community = best_communities[i]
		gt_comms_cnt = 0
		other_comms_cnt = 0
		for u in p:
			gt_comms = set(ground_truth.subsetsOf(u))
			if best_community not in gt_comms:
				if len(truth_communities.intersection(gt_comms)) > 0:
					gt_comms_cnt += 1
				else:
					other_comms_cnt += 1
		result_list.append({
			'partition_size': len(p),
			'wrong_nodes': gt_comms_cnt + other_comms_cnt,
			'wrong_nodes_gt': gt_comms_cnt,
			'wrong_nodes_other': other_comms_cnt,
			'wrong_percentage': float(gt_comms_cnt + other_comms_cnt) / num_partitions,
			'wrong_percentage_gt': float(gt_comms_cnt) / num_partitions,
			'wrong_percentage_other': float(other_comms_cnt) / num_partitions,
		})
		if gt_comms_cnt > 0.8 * len(p) and len(p) > 20:
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
		for i in range(len(partitions)):
			if i != best_partition:
				for u in partitions[i]:
					comms = ground_truth.subsetsOf(u)
					if c in comms:
						if best_communities[i] not in comms or best_communities[i] == c:
							node_cnt += 1

		result_list.append({
			'comm_size': community_sizes[c],
			'wrong_nodes': node_cnt,
			'wrong_percentage': float(node_cnt) / community_sizes[c],
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
	for i in range(len(partitions)):
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
