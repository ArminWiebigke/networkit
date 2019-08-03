import os
from collections import defaultdict
from copy import copy

from cover_benchmark import CoverBenchmark
from networkit import none
from networkit.components import ConnectedComponents
from egosplit.benchmarks.evaluation.output import *
from networkit.stopwatch import clockit


# Evaluate the partition of the ego-nets
def analyse_ego_net_partitions(benchmarks, result_dir, append, write_scores_per_egonet):
	open_mode = 'w'
	if append:
		open_mode = 'a'
	out_file_comm = open(result_dir + 'ego_net_communities.result', open_mode)
	out_file_part = open(result_dir + 'ego_net_partitions.result', open_mode)
	if write_scores_per_egonet:
		out_file_ego_metrics = open(result_dir + 'ego_net_ego_metrics.result', open_mode)
	else:
		out_file_ego_metrics = open(os.devnull, 'w')
	out_file_metrics = open(result_dir + 'ego_net_metrics.result', open_mode)
	out_file_conduct = open(result_dir + 'ego_net_conductance.result', open_mode)
	if not append:
		print_headers(out_file_comm, out_file_part, out_file_ego_metrics,
		              out_file_metrics, out_file_conduct)

	for benchmark in benchmarks:
		analyse_ego_net_partition(benchmark, out_file_comm, out_file_part,
		                          out_file_ego_metrics,
		                          out_file_metrics, out_file_conduct)


# Print the output file headers
def print_headers(out_file_comm, out_file_part, out_file_ego_metrics, out_file_metrics,
                  out_file_conduct):
	out_file_comm.write(create_line(
		*CoverBenchmark.output_header(),
		"comm_size",
		"wrong_nodes",
		"wrong_percentage",
		"num_partitions",
	))
	out_file_part.write(create_line(
		*CoverBenchmark.output_header(),
		"Ego-Node",
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
		*CoverBenchmark.output_header(),
		"Ego-Node",
		"Ego-Net Size",
		"Number of GT Communities",
		"Metric",
		"Value",
	))
	out_file_metrics.write(create_line(
		*CoverBenchmark.output_header(),
		"Metric",
		"Value",
	))
	out_file_conduct.write(create_line(
		*CoverBenchmark.output_header(),
		"Ego-Node",
		"Community Size",
		"Cut",
		"Volume",
		"Conductance"
	))


# Analyse the result of one benchmark run
def analyse_ego_net_partition(benchmark, out_comm, out_part, out_ego_metrics,
                              out_metrics, out_conduct):
	try:
		benchmark.algo.ego_net_partition_of(0)
	except AttributeError:
		return
	ground_truth = benchmark.get_ground_truth()
	graph = benchmark.get_graph()
	total_sums = defaultdict(lambda: 0)
	avg_metrics = defaultdict(lambda: 0)
	gt_comm_sizes = ground_truth.subsetSizeMap()
	for u in graph.nodes():
		truth_communities = set(ground_truth.subsetsOf(u))
		ego_net_size = graph.degree(u)
		ego_sums = defaultdict(lambda: 0)
		ego_net = benchmark.algo.ego_net_of(u)
		if ego_net.numberOfNodes() == 0:
			continue

		# Count nodes
		node_partition_map = copy(benchmark.algo.ego_net_partition_of(u))
		nodes = list(node_partition_map.keys())
		ego_sums["extended_ego_net_size"] = len(nodes)
		ego_sums["external_nodes_extended"] = count_external_nodes(
			ground_truth, truth_communities, nodes)

		# Evaluate structure of extended ego-net
		result_egonet_struct = eval_egonet_structure(
			ground_truth, truth_communities, ego_net)
		ego_sums["ext_intra_edges"] = result_egonet_struct["intra_edges"]
		ego_sums["ext_inter_edges"] = result_egonet_struct["inter_edges"]
		ego_sums["ext_egonet_edges"] = result_egonet_struct["egonet_edges"]
		# Per community
		result_list_comms_struct = eval_comms_structure(
			ground_truth, truth_communities, ego_net)
		ego_sums["ext_num_gt_comms"] = len(result_list_comms_struct)
		for result in result_list_comms_struct:
			ego_sums["ext_conductance"] += result["conductance"]
			ego_sums["ext_comm_fitness_0.5"] += result["comm_fitness_0.5"]
			ego_sums["ext_comm_fitness_0.6"] += result["comm_fitness_0.6"]
			ego_sums["ext_comm_fitness_0.7"] += result["comm_fitness_0.7"]
			ego_sums["ext_comm_fitness_0.8"] += result["comm_fitness_0.8"]
			ego_sums["ext_cut"] += result["cut"]
			ego_sums["ext_cut_comm"] += result["cut_comm"]
			ego_sums["ext_volume"] += result["volume"]
		# out_conduct.write(create_line(
		# 	algo_name,
		# 	graph_name,
		# 	u,
		# 	result["comm_size"],
		# 	result["cut"],
		# 	result["volume"],
		# 	result["conductance"],
		# ))

		# TODO: Modularity für Partition: GT comms, alle externen in einzelner comm

		# Number of connected compontents, nodes not in the largest connected component
		result_list_comm_comp = comm_components(ground_truth, truth_communities, ego_net)
		for result in result_list_comm_comp:
			ego_sums["num_components"] += result["num_components"]
			ego_sums["separate_nodes"] += result["separate_nodes"]

		# Coverage of the ground-truth communities
		result_coverage = calc_coverage(ground_truth, truth_communities, ego_net,
		                                gt_comm_sizes)
		for result in result_coverage:
			ego_sums["coverage"] += result["coverage"]

		# ********************************************************************************
		# *                Get original ego-net (remove extended nodes)                  *
		# ********************************************************************************
		num_extended_nodes = 0
		for v in nodes:
			if not graph.hasEdge(u, v):
				num_extended_nodes += 1
				del node_partition_map[v]
				ego_net.removeNode(v)
		nodes = list(node_partition_map.keys())
		ego_sums["extended_nodes"] = num_extended_nodes
		partitions, best_communities, community_sizes = calculate_partition_properties(
			ground_truth, truth_communities, node_partition_map)
		truth_communities = set([c for c in truth_communities if community_sizes[c] > 0])
		if not truth_communities:
			continue
		num_communities = len(truth_communities)
		ego_sums["num_comms"] = num_communities
		ego_sums["external_nodes"] = count_external_nodes(
			ground_truth, truth_communities, nodes)
		ego_sums["ego_net_size"] = ego_net_size

		# Evaluate structure of original ego-net
		result_egonet_struct = eval_egonet_structure(
			ground_truth, truth_communities, ego_net)
		ego_sums["intra_edges"] = result_egonet_struct["intra_edges"]
		ego_sums["inter_edges"] = result_egonet_struct["inter_edges"]
		ego_sums["egonet_edges"] = result_egonet_struct["egonet_edges"]
		# Per community
		result_list_comms_struct = eval_comms_structure(
			ground_truth, truth_communities, ego_net)
		ego_sums["num_gt_comms"] = len(result_list_comms_struct)
		for result in result_list_comms_struct:
			ego_sums["conductance"] += result["conductance"]
			ego_sums["cut"] += result["cut"]
			ego_sums["cut_comm"] += result["cut_comm"]
			ego_sums["volume"] += result["volume"]
		# out_conduct.write(create_line(
		# 	algo_name,
		# 	graph_name,
		# 	u,
		# 	result["comm_size"],
		# 	result["cut"],
		# 	result["volume"],
		# 	result["conductance"],
		# ))

		# Are communities split into multiple partitions?
		# TODO: Anzahl Communities, die eine Partition dominieren -> sinnvolle Personas
		# TODO: +1 wenn p beste Part. von Com. c und c dominiert p (Optimalfall)
		result_list_persona_recall = check_persona_recall(
			ground_truth, partitions, best_communities, community_sizes, truth_communities)
		for result in result_list_persona_recall:
			ego_sums["persona_recall"] += result["persona_recall"]

		result_list_comm = check_community_partition(
			ground_truth, partitions, community_sizes, truth_communities)
		for result in result_list_comm:
			ego_sums["comm_size"] += result["comm_size"]
			ego_sums["comm_incorrect"] += result["wrong_nodes"]
			ego_sums["parts_per_comm"] += result["num_partitions"]
		# out_comm.write(create_line(
		# 	algo_name,
		# 	graph_name,
		# 	result["comm_size"],
		# 	result["wrong_nodes"],
		# 	result["wrong_percentage"],
		# 	result["num_partitions"],
		# ))

		# Do partitions consist of multiple communities?
		result_list_part = check_partition_composition(
			ground_truth, partitions, best_communities, truth_communities)
		ego_sums["num_partitions"] = len(result_list_part)
		for result in result_list_part:
			ego_sums["partition_size"] += result["partition_size"]
			ego_sums["part_incorrect_gt"] += result["wrong_nodes_gt"]
			ego_sums["part_incorrect_ext"] += result["wrong_nodes_external"]
			ego_sums["comms_per_part"] += result["num_communities"]
		# out_part.write(create_line(
		# 	algo_name,
		# 	graph_name,
		# 	u,
		# 	result["partition_size"],
		# 	result["wrong_nodes"],
		# 	result["wrong_nodes_gt"],
		# 	result["wrong_nodes_external"],
		# 	result["wrong_percentage"],
		# 	result["wrong_percentage_gt"],
		# 	result["wrong_percentage_external"],
		# 	result["num_communities"],
		# ))

		# Count number of partitions
		result = count_partitions(partitions, truth_communities)
		for key, value in result.items():
			out_ego_metrics.write(create_line(
				*benchmark.output_line(),
				u,
				ego_net_size,
				num_communities,
				key,
				value,
			))

		# Summary metrics
		# TODO: Calculate scores per community/partition, take the average?
		#   makes small communities more/equally important
		for key, value in ego_sums.items():
			total_sums[key] += value
		ego_metrics = calc_metrics(ego_sums)
		for key, value in ego_metrics.items():
			avg_metrics[key] += value
		avg_metrics['num_values'] += 1

		for key, value in ego_metrics.items():
			out_ego_metrics.write(create_line(
				*benchmark.output_line(),
				u,
				ego_net_size,
				num_communities,
				key,
				value,
			))

	total_metrics = {key: (value / avg_metrics['num_values'])
	                 for (key, value) in avg_metrics.items()}
	# total_metrics = calc_metrics(total_sums)

	for key, value in total_metrics.items():
		out_metrics.write(create_line(
			*benchmark.output_line(),
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
		"Community Segmentation": sums["comm_incorrect"] / sums["comm_size"],
		"Cluster per Community": sums["parts_per_comm"] / sums["num_comms"],
		"Communities per Cluster": safe_div(sums["comms_per_part"], sums["num_partitions"], 1),
		"Merged External Nodes": safe_div(sums["part_incorrect_ext"],
		                                  sums["partition_size"]),
		"Community Merging": safe_div(sums["part_incorrect_gt"],
		                                  sums["partition_size"]),
		"Extended Nodes / Ego-Net Size": safe_div(sums["extended_nodes"], sums["ego_net_size"]),
		"External Nodes Ratio": safe_div(sums["external_nodes_extended"],
		                           sums["extended_ego_net_size"]),
		"Added External Nodes Ratio": safe_div(
			sums["external_nodes_extended"] - sums["external_nodes"],
			sums["extended_ego_net_size"] - sums["ego_net_size"]),
		"Conductance": safe_div(sums["ext_conductance"], sums["ext_num_gt_comms"]),
		# "conductance_ratio": safe_div(
		# 	(sums["ext_conductance"] / sums["ext_num_gt_comms"]),
		# 	(sums["conductance"] / sums["num_gt_comms"]), 1),
		"Community Fitness": safe_div(sums["ext_comm_fitness_0.8"], sums["ext_num_gt_comms"]),
		"Ratio of Intra-Community Edges": safe_div(sums["ext_intra_edges"], sums["ext_egonet_edges"]),
		"Change of Intra-Edges": safe_div(
			safe_div(sums["ext_intra_edges"], sums["ext_egonet_edges"]),
			safe_div(sums["intra_edges"], sums["egonet_edges"]), 1) - 1,
		"Components per Community": safe_div(sums["num_components"], sums["ext_num_gt_comms"]),
		"Ratio of Disconnected Nodes": safe_div(sums["separate_nodes"], sums["ext_num_gt_comms"]),
		"Persona Recall": safe_div(sums["persona_recall"], sums["num_gt_comms"]),
		"Community Coverage": safe_div(sums["coverage"], sums["num_gt_comms"]),
		"Number of Personas": sums["num_partitions"],
		'Number of Ground-Truth Communities': sums["num_gt_comms"],
	}

	# a = 1 - metrics["Community Segmentation"]
	# b = 1 - metrics["Community Merging"]
	# metrics["ego_partition_score_harm"] = 1 - harmonic_mean(a, b)
	# metrics["ego_partition_score_arit"] = (metrics["community_cohesion"] +
	#                                        metrics["partition_exclusivity"]) / 2
	return metrics


def calc_coverage(ground_truth, truth_communities, ego_net, gt_comm_sizes):
	gt_comms, _ = get_gt_communities(ego_net, ground_truth, truth_communities)
	result_list = []
	for c in truth_communities:
		coverage = len(gt_comms[c]) / gt_comm_sizes[c]
		result_list.append({
			'coverage': coverage,
		})
	return result_list


def harmonic_mean(a, b):
	if a + b == 0:
		return 0
	return 2 * a * b / (a + b)


def check_persona_recall(ground_truth, partitions, best_communities, community_sizes,
                         truth_communities):
	result_list = []
	best_partitions = {}
	best_part_size = {}

	for c in truth_communities:
		best, s = find_best_partition(ground_truth, c, partitions)
		best_partitions[c] = best
		best_part_size[c] = s
	for c in truth_communities:
		best_node_cnt = 0
		for i, partition in enumerate(partitions):
			if best_communities[i] == c:
				node_cnt = count_comm_nodes_in_partition(c, ground_truth, partition)
				best_node_cnt = max(best_node_cnt, node_cnt)
		persona_recall = best_node_cnt / community_sizes[c]

		result_list.append({
			'comm_size': community_sizes[c],
			'persona_recall': persona_recall,
		})
	return result_list


def count_comm_nodes_in_partition(c, ground_truth, partition):
	num_nodes = 0
	for u in partition:
		if c in ground_truth.subsetsOf(u):
			num_nodes += 1
	return num_nodes


# Count the number of connected components for each ground truth community
# @clockit
def comm_components(ground_truth, truth_communities, ego_net):
	communities, _ = get_gt_communities(ego_net, ground_truth, truth_communities)

	result_list = []
	for comm, nodes in communities.items():
		subgraph = ego_net.subgraphFromNodes(nodes)
		algo = ConnectedComponents(subgraph)
		algo.run()
		num_components = algo.numberOfComponents()
		component_sizes = algo.getComponentSizes()
		max_component, max_size = max(component_sizes.items())

		result_list.append({
			'comm_size': len(nodes),
			'num_components': num_components,
			'separate_nodes': (len(nodes) - max_size) / len(nodes)
		})
	return result_list


# TODO: Auch externe Knoten können die Community stärker verbinden
# TODO: Warum wird das Ergebnis mit vielen externen Knoten besser? (edges vs. sign)

# Count the intra/inter community edges on the whole ego-net
# @clockit
def eval_egonet_structure(ground_truth, truth_communities, ego_net):
	num_intra = 0
	for u, v in ego_net.edges():
		common = ground_truth.subsetsOf(u).intersection(ground_truth.subsetsOf(v))
		if common.intersection(truth_communities):
			num_intra += 1
	num_inter = ego_net.numberOfEdges() - num_intra
	return {
		'intra_edges': num_intra,
		'inter_edges': num_inter,
		'egonet_edges': num_intra + num_inter,
	}


# TODO: Add external nodes to commuities (if at least 50% of edges to one community),
#   then evaluate community structure
# Calculate conductance of each community
# @clockit
def eval_comms_structure(ground_truth, truth_communities, ego_net):
	communities, external_nodes = get_gt_communities(ego_net, ground_truth,
	                                                 truth_communities)
	external_nodes = set(external_nodes)

	result_list = []
	for comm, nodes in communities.items():
		cut_comm = 0
		cut_ext = 0
		volume = 0
		for u in nodes:
			volume += ego_net.degree(u)
			for v in ego_net.neighbors(u):
				if not v in nodes:
					if v in external_nodes:
						cut_ext += 1
					else:
						cut_comm += 1

		cut = cut_comm + cut_ext
		inner_degree = volume - cut
		if cut == 0:
			conductance = 0
		else:
			conductance = cut / volume
		if volume == 0:
			fitness_0_5, fitness_0_6, fitness_0_7, fitness_0_8 = 0, 0, 0, 0
		else:
			fitness_0_5 = inner_degree / (volume ** 0.5)
			fitness_0_6 = inner_degree / (volume ** 0.6)
			fitness_0_7 = inner_degree / (volume ** 0.7)
			fitness_0_8 = inner_degree / (volume ** 0.8)

		result_list.append({
			'comm_size': len(nodes),
			'cut': cut,
			'cut_comm': cut_comm,
			'volume': volume,
			'conductance': conductance,
			'comm_fitness_0.5': fitness_0_5,
			'comm_fitness_0.6': fitness_0_6,
			'comm_fitness_0.7': fitness_0_7,
			'comm_fitness_0.8': fitness_0_8,
		})
	return result_list


# Return a map that maps the ground truth communities to its nodes in the ego-net
# gt_id -> list of nodes
def get_gt_communities(ego_net, ground_truth, truth_communities):
	communities = defaultdict(lambda: [])
	external_nodes = []
	for u in ego_net.nodes():
		comms = ground_truth.subsetsOf(u).intersection(truth_communities)
		if not comms:
			external_nodes.append(u)
		for c in comms:
			communities[c].append(u)
	return communities, external_nodes


# Count the number of nodes that don't share any communities with the ego-node
def count_external_nodes(ground_truth, truth_communities, nodes):
	cnt = 0
	for v in nodes:
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
	community_cnts = []  # for each partition: number of nodes for each ground truth community
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
# @clockit
def check_partition_composition(ground_truth, partitions, best_communities,
                                truth_communities):
	result_list = []
	for p_id, partition in enumerate(partitions):
		partition_size = len(partition)
		best_community = best_communities[p_id]
		if best_community == -1:
			continue
		gt_comms_cnt = 0
		external_nodes = 0
		community_nodes = defaultdict(lambda: 0)
		for u in partition:
			gt_comms = set(ground_truth.subsetsOf(u))
			if best_community not in gt_comms:
				# Found a node that is not in the dominating community
				comms = truth_communities.intersection(gt_comms)
				if len(comms) > 0:
					gt_comms_cnt += 1
					for c in comms:
						community_nodes[c] += 1
				else:
					external_nodes += 1
		# Count overlapping nodes only once, for their largest community
		for u in partition:
			gt_comms = set(ground_truth.subsetsOf(u)).intersection(truth_communities)
			if not gt_comms or best_community in gt_comms:
				continue
			comm_counts = [(community_nodes[c], c) for c in gt_comms]
			_, max_comm = max(comm_counts)
			for c in gt_comms:
				if c != max_comm:
					community_nodes[c] -= 1
		num_communities = 1 + len([c for c in community_nodes.values() if c > 0])

		result_list.append({
			'partition_size': partition_size,
			'wrong_nodes': gt_comms_cnt + external_nodes,
			'wrong_nodes_gt': gt_comms_cnt,
			'wrong_nodes_external': external_nodes,
			'wrong_percentage': 100 * (
					gt_comms_cnt + external_nodes) / partition_size,
			'wrong_percentage_gt': 100 * gt_comms_cnt / partition_size,
			'wrong_percentage_external': 100 * external_nodes / partition_size,
			'num_communities': num_communities,
		})
	return result_list


# Count the number of nodes of the given community in all partitions except the partition
# with the largest number of nodes from that community.
# Don't count nodes that are in multiple communities and that are assigned to one of its
# communities correctly.
# @clockit
def check_community_partition(ground_truth, partitions, community_sizes,
                              truth_communities):
	result_list = []
	best_partitions = {}
	for c in truth_communities:
		best_partitions[c], _ = find_best_partition(ground_truth, c, partitions)
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
		num_nodes = count_comm_nodes_in_partition(community, ground_truth, partition)
		if num_nodes > best_partition_cnt:
			best_partition_cnt = num_nodes
			best_partition = i
	return best_partition, best_partition_cnt


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
