from collections import defaultdict
from copy import copy

from networkit import none
from egosplit.benchmarks.evaluation.output import *


# Evaluate the partition of the ego-nets
def analyse_ego_net_partitions(benchmarks, result_dir, append):
	open_mode = 'w'
	if append:
		open_mode = 'a'
	out_file_comm = open(result_dir + 'ego_net_communities.result', open_mode)
	out_file_part = open(result_dir + 'ego_net_partitions.result', open_mode)
	out_file_ego_metrics = open(result_dir + 'ego_net_ego_metrics.result', open_mode)
	out_file_metrics = open(result_dir + 'ego_net_metrics.result', open_mode)
	if not append:
		print_headers(out_file_comm, out_file_part, out_file_ego_metrics, out_file_metrics)

	for benchmark in benchmarks:
		analyse_ego_net_partition(benchmark, out_file_comm, out_file_part, out_file_ego_metrics,
		                          out_file_metrics)


# Print the output file headers
def print_headers(out_file_comm, out_file_part, out_file_ego_metrics, out_file_metrics):
	out_file_comm.write(create_line(
		"algo",
		"graph",
		"comm_size",
		"wrong_nodes",
		"wrong_percentage",
		"num_partitions",
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
		"wrong_percentage_other",
		"num_communities",
	))
	out_file_ego_metrics.write(create_line(
		"algo",
		"graph",
		"node",
		"ego_net_size",
		"num_comms",
		"metric_name",
		"value",
	))
	out_file_metrics.write(create_line(
		"algo",
		"graph",
		"metric_name",
		"value",
	))


# Analyse the result of one benchmark run
def analyse_ego_net_partition(benchmark, out_comm, out_part, out_ego_metrics,
                              out_metrics):
	try:
		ego_net_partitions = benchmark.algo.ego_net_partition_of(0)
	except AttributeError:
		return
	ground_truth = benchmark.get_ground_truth()
	algo_name = benchmark.algo.name
	graph_name = benchmark.graph.name
	graph = benchmark.get_graph()
	total_sums = defaultdict(lambda: 0)
	for u in graph.nodes():
		truth_communities = set(ground_truth.subsetsOf(u))
		ego_net_size = graph.degree(u)
		ego_sums = defaultdict(lambda: 0)

		# Remove all nodes that are not neighbors of the ego-node
		node_partition_map = copy(benchmark.algo.ego_net_partition_of(u))
		nodes = list(node_partition_map.keys())
		ego_sums["extended_ego_net_size"] = len(nodes)
		ego_sums["external_nodes_extended"] = count_external_nodes(
			ground_truth, truth_communities, nodes)
		num_extended_nodes = 0
		for v in nodes:
			if not graph.hasEdge(u, v):
				num_extended_nodes += 1
				del node_partition_map[v]

		partitions, best_communities, community_sizes = calculate_partition_properties(
			ground_truth, truth_communities, node_partition_map)
		truth_communities = set([c for c in truth_communities if community_sizes[c] > 0])
		if not truth_communities:
			continue
		num_communities = len(truth_communities)
		ego_sums["extended_nodes"] = num_extended_nodes
		ego_sums["external_nodes"] = count_external_nodes(
			ground_truth, truth_communities, node_partition_map.keys())
		ego_sums["ego_net_size"] = ego_net_size

		# Are communities split into multiple partitions?
		# TODO: Anzahl Communities, die eine Partition dominieren -> sinnvolle Personas
		# TODO: +1 wenn p beste Part. von Com. c und c dominiert p (Optimalfall)
		result_list_comm = check_community_partition(
			ground_truth, partitions, community_sizes, truth_communities)
		ego_sums["num_comms"] = len(truth_communities)
		for result in result_list_comm:
			ego_sums["comm_size"] += result["comm_size"]
			ego_sums["comm_incorrect"] += result["wrong_nodes"]
			ego_sums["parts_per_comm"] += result["num_partitions"]
			out_comm.write(create_line(
				algo_name,
				graph_name,
				result["comm_size"],
				result["wrong_nodes"],
				result["wrong_percentage"],
				result["num_partitions"],
			))

		# Do partitions consist of multiple communities?
		result_list_part = check_partition_composition(
			ground_truth, partitions, best_communities, truth_communities)
		ego_sums["num_partitions"] = len(result_list_part)
		for result in result_list_part:
			ego_sums["partition_size"] += result["partition_size"]
			ego_sums["part_incorrect_gt"] += result["wrong_nodes_gt"]
			ego_sums["part_incorrect_ext"] += result["wrong_nodes_external"]
			ego_sums["comms_per_part"] += result["num_communities"]
			out_part.write(create_line(
				algo_name,
				graph_name,
				result["partition_size"],
				result["wrong_nodes"],
				result["wrong_nodes_gt"],
				result["wrong_nodes_external"],
				result["wrong_percentage"],
				result["wrong_percentage_gt"],
				result["wrong_percentage_external"],
				result["num_communities"],
			))

		# Count number of partitions
		result = count_partitions(partitions, truth_communities)
		for key, value in result.items():
			out_ego_metrics.write(create_line(
				algo_name,
				graph_name,
				u,
				ego_net_size,
				num_communities,
				key,
				value,
			))

		# Summary metrics
		ego_metrics = calc_metrics(ego_sums)
		for key, value in ego_metrics.items():
			out_ego_metrics.write(create_line(
				algo_name,
				graph_name,
				u,
				ego_net_size,
				num_communities,
				key,
				value,
			))
		for key, value in ego_sums.items():
			total_sums[key] += value

		# TODO: Count external nodes -> Good extension?

	total_metrics = calc_metrics(total_sums)

	for key, value in total_metrics.items():
		out_metrics.write(create_line(
			algo_name,
			graph_name,
			key,
			value,
		))


# Calculate the ego-net metrics
def calc_metrics(sums):
	def safe_div(a, b, default=0):
		if b == 0:
			return default
		return a / b
	metrics = {
		"community_cohesion": sums["comm_incorrect"] / sums["comm_size"],
		"parts_per_comm": sums["parts_per_comm"] / sums["num_comms"],
		"comms_per_part": safe_div(sums["comms_per_part"], sums["num_partitions"], 1),
		"merged_external_nodes": safe_div(sums["part_incorrect_ext"],
		                                             sums["external_nodes"]),
		"partition_exclusivity": safe_div(sums["part_incorrect_gt"],
		                                              sums["partition_size"]),
		"extended_nodes": safe_div(sums["extended_nodes"], sums["ego_net_size"]),
		"external_nodes": safe_div(sums["external_nodes_extended"],
		                                    sums["extended_ego_net_size"]),
	}
	a = 1 - metrics["community_cohesion"]
	b = 1 - metrics["partition_exclusivity"]
	metrics["ego_partition_score"] = 1 - 2 * a * b / (a + b)
	return metrics


# Count the number of nodes that don't share any communities with the ego-node
def count_external_nodes(ground_truth, truth_communities, neighbors):
	cnt = 0
	for v in neighbors:
		if not truth_communities.intersection(ground_truth.subsetsOf(v)):
			cnt += 1
	return cnt


# Returns the number of partitions in the ego-net
def count_partitions(partitions, truth_communitites):
	num_partitions = len(partitions)
	num_two_plus_partitions = len([p for p in partitions if len(p) >= 2])
	num_three_plus_partitions = len([p for p in partitions if len(p) >= 3])
	num_truth_communities = len(truth_communitites)

	result = {
		'num_partitions': num_partitions,
		'num_two_plus_partitions': num_two_plus_partitions,
		'num_three_plus_partitions': num_three_plus_partitions,
		'truth_communities': num_truth_communities,
	}

	return result


# Returns partition->node list, best community for each partition, community->size map
def calculate_partition_properties(ground_truth, truth_communities, partition_map):
	max_partition = 0
	for id in partition_map.values():
		max_partition = max(max_partition, id)
	partitions = [[] for _ in range(max_partition + 1)]
	for v, p_id in partition_map.items():
		if v != none and p_id != none:  # TODO: Why does leidenalg assign partition = none, but no other algorithm?
			partitions[p_id].append(v)
	# Remove empty partitions
	partitions = [p for p in partitions if len(p) > 0]
	max_community_cnts = {}
	for c in truth_communities:
		max_community_cnts[c] = 0
	community_sizes = defaultdict(lambda: 0)
	best_communities = []  # for each partition: the ground truth community with the most nodes
	community_cnts = []    # for each partition: number of nodes for each ground truth community
	for p_id, partition in enumerate(partitions):
		community_cnts.append(defaultdict(lambda: 0))
		for v in partition:
			comms = truth_communities.intersection(ground_truth.subsetsOf(v))
			for c in comms:
				community_cnts[p_id][c] += 1
		for c in community_cnts[p_id]:
			max_community_cnts[c] = max(max_community_cnts[c], community_cnts[p_id][c])
			community_sizes[c] += community_cnts[p_id][c]

		best_communities.append(dominating_community(community_cnts[p_id]))
	return partitions, best_communities, community_sizes


# For each partition: Get the community with the most nodes. Count the number of nodes
# that are not in that community.
def check_partition_composition(ground_truth, partitions, best_communities,
                                truth_communities):
	result_list = []
	for p_id, partition in enumerate(partitions):
		partition_size = len(partition)
		best_community = best_communities[p_id]
		if best_community == -1:
			continue
		gt_comms_cnt = 0
		external_comms_cnt = 0
		community_in_partition = {}
		for u in partition:
			gt_comms = set(ground_truth.subsetsOf(u))
			if best_community not in gt_comms:
				# Found a node that is not in the dominating community
				comms = truth_communities.intersection(gt_comms)
				if len(comms) > 0:
					gt_comms_cnt += 1
					for c in comms:
						community_in_partition[c] = True
				else:
					external_comms_cnt += 1
		num_communities = 1 + len(community_in_partition)  # TODO: Overlapping nodes count multiple times
		result_list.append({
			'partition_size': partition_size,
			'wrong_nodes': gt_comms_cnt + external_comms_cnt,
			'wrong_nodes_gt': gt_comms_cnt,
			'wrong_nodes_external': external_comms_cnt,
			'wrong_percentage': 100 * (gt_comms_cnt + external_comms_cnt) / partition_size,
			'wrong_percentage_gt': 100 * gt_comms_cnt / partition_size,
			'wrong_percentage_external': 100 * external_comms_cnt / partition_size,
			'num_communities': num_communities,
		})
		return result_list
	return result_list


# Count the number of nodes of the given community in all partitions except the partition
# with the largest number of nodes from that community.
# Don't count nodes that are in multiple communities and that are assigned to one of its
# communities correctly.
def check_community_partition(ground_truth, partitions, community_sizes, truth_communities):
	result_list = []
	best_partitions = {}
	for c in truth_communities:
		best_partitions[c] = find_best_partition(ground_truth, c, partitions)
	for c in truth_communities:
		best_partition = best_partitions[c]
		node_cnt = 0
		node_in_partition = [0 for _, _ in enumerate(partitions)]
		for p_id, partition in enumerate(partitions):
			if p_id == best_partition:
				continue
			for u in partition:
				comms = truth_communities.intersection(ground_truth.subsetsOf(u))
				if c in comms:
					# Found a node of the community in another partition.
					# Don't count the node if it has another ground truth community that has most
					# of it's nodes in this partition (relative majority).
					count_node = True
					for other_c in comms:
						if other_c != c and best_partitions[other_c] == p_id:
							count_node = False
					if count_node:
						node_cnt += 1
						node_in_partition[p_id] = 1
		num_partitions = 1 + len([x for x in node_in_partition if x > 0])

		result_list.append({
			'comm_size': community_sizes[c],
			'wrong_nodes': node_cnt,
			'wrong_percentage': 100 * node_cnt / community_sizes[c],
			'num_partitions': num_partitions
		})
	return result_list


# Find the partition with the largest number of nodes that are in the given community
def find_best_partition(ground_truth, community, partitions):
	best_partition = -1
	best_partition_cnt = 0
	for i, partition in enumerate(partitions):
		num_nodes = 0
		for u in partition:
			if community in ground_truth.subsetsOf(u):
				num_nodes += 1
		if num_nodes > best_partition_cnt:
			best_partition_cnt = num_nodes
			best_partition = i
	return best_partition


# For a given partition, return the community with the highest number of nodes.
# If no community has two or more nodes, return -1.
def dominating_community(community_cnts):
	max_comm = -1
	max_comm_cnt = 0
	for comm, comm_cnt in community_cnts.items():
		if comm_cnt > max_comm_cnt and comm_cnt >= 2:
			max_comm = comm
			max_comm_cnt = comm_cnt
	return max_comm
