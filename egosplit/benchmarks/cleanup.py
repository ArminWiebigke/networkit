from copy import copy
from random import randint
from collections import OrderedDict, defaultdict

from networkit.community import CoverF1Similarity
from networkit.structures import Cover
from networkit import none
from egosplit.external import calc_entropy


def poss_combs(x):
	return x * (x - 1) / 2


# Merge highly connected partitions
def partition_merge(G, input_partition):
	# return input_partition
	result = copy(input_partition)
	continue_merge = True
	# print(G.numberOfNodes())
	# print(poss_combs(G.numberOfNodes()))
	graph_density = G.numberOfEdges() / poss_combs(G.numberOfNodes())
	# print(graph_density)  # TODO: This should be the density of the global graph instead
	while continue_merge:
		continue_merge = False
		partition_sets = [set() for _ in range(result.upperBound())]
		for u in G.nodes():
			part_id = result.subsetOf(u)
			partition_sets[part_id].add(u)
		partitions = [(i, p) for (i, p) in enumerate(partition_sets) if len(p) > 1]

		cuts = [defaultdict(lambda: 0) for _ in range(result.upperBound())]
		for u in G.nodes():
			part_id = result.subsetOf(u)
			for v in G.neighbors(u):
				cuts[part_id][result.subsetOf(v)] += 1

		partition_values = [defaultdict(lambda: 0.0) for _ in range(result.upperBound())]
		for part_id, partition in partitions:
			# print(part_id, partition)
			internal_degree = cuts[part_id][part_id] / 2
			part_size = len(partition)
			if part_size > 1:
				internal_density = internal_degree / poss_combs(part_size)
			else:
				internal_density = 1
			partition_values[part_id]["density"] = internal_density
			partition_values[part_id]["edges"] = internal_degree
			partition_values[part_id]["size"] = part_size

		for id_a, part_a in partitions:
			# print(id_a)
			# print(partition_values[id_a])
			density = partition_values[id_a]["density"]
			for id_b, part_b in partitions:
				# print("\t", id_b)
				if id_a == id_b:
					continue
				if cuts[id_a][id_b] <= 2:
					continue

				size_a = partition_values[id_a]["size"]
				size_b = partition_values[id_b]["size"]
				combined_size = size_a + size_b
				edges_a = partition_values[id_a]["edges"]
				edges_b = partition_values[id_b]["edges"]
				cut = cuts[id_a][id_b]
				combined_edges = edges_a + edges_b + cut
				merge = False

				cut_density = cut / (size_a * size_b)
				separate_density = (edges_a + edges_b) / (poss_combs(size_a) + poss_combs(size_b))
				combined_density = combined_edges / poss_combs(combined_size)
				merge = cut_density > 0.1
				# if cut > edges_a or cut > edges_b:
				# 	if combined_density > 1.5 * graph_density:
				# 		merge = True
				# if cut_density > graph_density:
				# 	merge = True
				# merge = cut / combined_edges > 0.2

				if merge:
					print("\t\t{}-{}: {}".format(size_a, size_b, cut), cut_density, separate_density, combined_density, graph_density, sep=", ")
					# print("############### Merge {} and {} ################".format(id_a, id_b))
					for u in part_b:
						result.moveToSubset(id_a, u)
					continue_merge = True
				if continue_merge:
					break
			if continue_merge:
				break
	return result


# Remove node from a community if this decreases the entropy.
def entropy_trim(G, cover):
	cover = copy(cover)
	communities = get_community_vector(G, cover)

	entropy = calc_entropy(G, cover)
	for c_id, community in enumerate(communities):
		comm_size = len(community)
		if comm_size == 0:
			continue
		print("Community:", c_id)
		community_changed = True
		while community_changed:
			community_changed = False
			not_removed_cnt = 0
			node_candidates = worst_nodes_in_comm(G, community)
			for _, u in node_candidates:
				cover.removeFromSubset(c_id, u)
				new_entropy = calc_entropy(G, cover)
				if new_entropy < entropy:
					# Remove node
					print("Remove node ({})".format(new_entropy-entropy))
					entropy = new_entropy
					community.remove(u)
					comm_size -= 1
					community_changed = True
					break
				else:
					# Leave node in the community
					cover.addToSubset(c_id, u)
					not_removed_cnt += 1
					if not_removed_cnt / comm_size > 0.10 and not_removed_cnt >= 5:
						break

	return cover


# Remove communities if this decreases the entropy.
def remove_comms_entropy(G, cover):
	cover = copy(cover)
	communities = get_community_vector(G, cover)

	entropy = calc_entropy(G, cover)
	for c_id, community in enumerate(communities):
		print("Community:", c_id)
		if len(community) >= 50 or len(community) == 0:
			continue
		for u in community:
			cover.removeFromSubset(c_id, u)
		new_entropy = calc_entropy(G, cover)
		if new_entropy < entropy:
			print("Remove community {} of size {} ({})".format(c_id, len(community),
			                                                   new_entropy-entropy))
			entropy = new_entropy
		else:
			for u in community:
				cover.addToSubset(c_id, u)

	return cover


# TODO: Merge und Trim verbinden?
# Merge communities together if they are strongly connected.
# def merge_comms_entropy_from_cover(G, cover):
# 	cover = copy(cover)
# 	communities = get_community_vector(G, cover)
# 	comm_sizes = [len(c) for c in communities]
# 	size_map = cover.subsetSizeMap()
# 	for i, c in enumerate(communities):
# 		if len(c) == 0:
# 			continue
# 		assert(size_map[i] == comm_sizes[i])
# 		assert(size_map[i] == len(c))
# 	print(comm_sizes)
# 	entropy = calc_entropy(G, cover)
# 	for c_id, community in enumerate(communities):
# 		max_comm_size = 50
# 		if len(community) >= max_comm_size or len(community) == 0:
#           # Only for the specific LFR graph
# 			continue
# 		community_changed = True
# 		while community_changed:
# 			community_changed = False
# 			cuts = community_cut(G, cover, comm_sizes, community, c_id)
# 			iterations = 3
# 			for cut_score, other_id in cuts:
# 				assert(cut_score <= 2.0)
# 				if cut_score < 0.07:  # Always useful?
# 					break
# 				if comm_sizes[other_id] > max_comm_size:
# 					continue
# 				if comm_sizes[other_id] < comm_sizes[c_id]:
# 					print("Skip {} -> {} ({})".format(
# 						comm_sizes[c_id], comm_sizes[other_id], cut_score))
# 					continue
# 				if iterations <= 0:
# 					break
# 				iterations -= 1
# 				print(cut_score)
# 				new_cover = copy(cover)
# 				other_comm = communities[other_id]
# 				for v in other_comm - community:
# 					new_cover.addToSubset(c_id, v)
# 				for v in other_comm:
# 					new_cover.removeFromSubset(other_id, v)
#
# 				new_entropy = calc_entropy(G, new_cover)
# 				overlap = len(communities[c_id].intersection(communities[other_id]))
# 				overlap /= min(comm_sizes[c_id], comm_sizes[other_id])
# 				if overlap >= 0.5:
# 					out_file = open("./results/comm_overlap.out", 'a')
# 					out_str = "IDs: {} {}, Sizes: {} {}, OL: {:.2f}, Score: {:.2f}, " \
# 					          "Entr:{:.2f}\n"
# 					out_file.write(out_str.format(
# 						c_id, other_id, comm_sizes[c_id], comm_sizes[other_id],
# 						overlap, cut_score, new_entropy - entropy
# 					))
# 					out_file.close()
# 				if new_entropy < entropy:
# 					print("\tMerge community with score {} #####".format(cut_score))
# 					entropy = new_entropy
# 					cover = new_cover
# 					comm_sizes[c_id] = len(community.union(other_comm))
# 					for u in other_comm:
# 						community.add(u)
# 					other_comm.clear()
# 					comm_sizes[other_id] = 0
# 					community_changed = True
# 					break
# 	return cover


# Merge communities together if they are strongly connected.
def merge_comms_entropy(G, cover, bad_groups):
	cover = copy(cover)
	good_comms = defaultdict(lambda: True)
	good_group_bound = cover.upperBound()
	bad_groups.sort(key=len)
	cover.setUpperBound(good_group_bound + len(bad_groups))
	bad_groups = list(enumerate(bad_groups, good_group_bound))

	# Add bad groups to cover
	for id, bad_group in bad_groups:
		good_comms[id] = False
		for u in bad_group:
			cover.addToSubset(id, u)

	communities = get_community_vector(G, cover)
	comm_sizes = [len(c) for c in communities]
	size_map = cover.subsetSizeMap()
	for i, c in enumerate(communities):
		if len(c) == 0:
			continue
		assert(size_map[i] == comm_sizes[i])
		assert(size_map[i] == len(c))
	print(comm_sizes)
	entropy = calc_entropy(G, cover)
	for c_id, community in enumerate(communities):
		if c_id < good_group_bound or len(community) == 0:
			continue
		community_changed = True
		while community_changed:
			community_changed = False
			cuts = community_cut(G, cover, comm_sizes, community, c_id)
			iterations = 3
			for cut_score, other_id in cuts:
				if cut_score < 0.1:  # Always useful?
					break
				if other_id < good_group_bound:
					continue
				if other_id >= c_id:
					print("Skip {} -> {} ({})".format(c_id, other_id, cut_score))
					continue
				# if comm_sizes[other_id] < comm_sizes[c_id]:
				# 	print("Skip {} -> {} ({})".format(
				# 		comm_sizes[c_id], comm_sizes[other_id], cut_score))
				# 	continue
				if iterations <= 0:
					break
				print("{} -> {} ({})".format(c_id, other_id, cut_score))
				iterations -= 1
				print(cut_score)
				new_cover = copy(cover)
				other_comm = communities[other_id]
				for v in other_comm - community:
					new_cover.addToSubset(c_id, v)
				for v in other_comm:
					new_cover.removeFromSubset(other_id, v)

				new_entropy = calc_entropy(G, new_cover)
				overlap = len(communities[c_id].intersection(communities[other_id]))
				overlap /= min(comm_sizes[c_id], comm_sizes[other_id])
				if overlap >= 0.5:
					out_file = open("./results/comm_overlap.out", 'a')
					out_str = "IDs: {} {}, Sizes: {} {}, OL: {:.2f}, Score: {:.2f}, " \
					          "Entr:{:.2f}\n"
					out_file.write(out_str.format(
						c_id, other_id, comm_sizes[c_id], comm_sizes[other_id],
						overlap, cut_score, new_entropy - entropy
					))
					out_file.close()
				if new_entropy < entropy:
					print("\tMerge community with score {} #####".format(cut_score))
					good_comms[c_id] = True
					good_comms[other_id] = True
					entropy = new_entropy
					cover = new_cover
					comm_sizes[c_id] = len(community.union(other_comm))
					for u in other_comm:
						community.add(u)
					other_comm.clear()
					comm_sizes[other_id] = 0
					community_changed = True
					break

	# Remove bad groups that were never merged
	for id, bad_group in bad_groups:
		if not good_comms[id]:
			for u in bad_group:
				cover.removeFromSubset(id, u)
	return cover


# For each found community, calculate the best fitting ground truth community.
# Then remove all nodes from the community that are not in the ground truth community.
def trim_comms(G, ground_truth, cover):
	cover = copy(cover)
	communities = get_community_vector(G, cover)
	gt_communities = get_community_vector(G, ground_truth)

	for c_id, community in enumerate(communities):
		_, best_gt_comm = find_best_gt_community(community, gt_communities)
		while remove_bad_node(G, cover, c_id, community, best_gt_comm):
			pass

	return cover


# Remove communities of size smaller than 5.
def remove_small_comms(G, cover, min_size=5):
	cover = copy(cover)

	communities = [set() for _ in range(cover.upperBound())]
	for u in G.nodes():
		for c in cover.subsetsOf(u):
			if len(communities[c]) < min_size:
				communities[c].add(u)

	for i, c in enumerate(communities):
		if len(c) < min_size:
			for u in c:
				cover.removeFromSubset(i, u)

	return cover


# Merge all communities that have the same ground truth community as their best fitting
# community.
def merge_gt_comms(G, ground_truth, cover):
	cover = copy(cover)
	communities = get_community_vector(G, cover)
	gt_communities = get_community_vector(G, ground_truth)

	comm_sets = defaultdict(lambda: [])  # The detected communities for each GT community
	for c_id, community in enumerate(communities):
		best_gt_comm, _ = find_best_gt_community(community, gt_communities)
		comm_sets[best_gt_comm].append(c_id)

	for merge_comms in comm_sets.values():
		# Merge all other communities into the first
		comm_a = merge_comms.pop(0)
		for comm_b in merge_comms:
			for u in communities[comm_b]:
				cover.addToSubset(comm_a, u)
				cover.removeFromSubset(comm_b, u)
	return cover


def remove_overlap_comms(G, cover, min_overlap=0.7):
	cover = copy(cover)
	communities = get_community_vector(G, cover)
	for comm_a, nodes_a in enumerate(communities):
		if len(nodes_a) == 0:
			continue
		for comm_b, nodes_b in enumerate(communities):
			if len(nodes_b) == 0 or comm_a == comm_b:
				continue
			common_nodes = nodes_a.intersection(nodes_b)
			overlap = len(common_nodes) / len(nodes_b)
			if overlap < min_overlap:
				continue
			for u in nodes_b:
				cover.removeFromSubset(comm_b, u)
			nodes_b.clear()
	return cover


# Merge a community A into another community B if most nodes of A are already in B
def merge_overlap_comms(G, cover, min_overlap=0.8):
	cover = copy(cover)
	communities = get_community_vector(G, cover)
	for a, comm_a in enumerate(communities):
		candidates = []
		for b, comm_b in enumerate(communities):
			if a == b:
				continue
			common_nodes = comm_a.intersection(comm_b)
			if not common_nodes:
				continue
			overlap = len(common_nodes) / min(len(comm_a), len(comm_b))
			if overlap >= min_overlap:
				candidates.append((overlap, b))
		# Merge with the community with the highest overlap
		if candidates:
			_, b = max(candidates)
			comm_b = communities[b]
			# Merge b into a and remove b
			for u in comm_b - comm_a:
				cover.addToSubset(a, u)
				comm_a.add(u)
			for u in comm_b:
				cover.removeFromSubset(b, u)
			comm_b.clear()
	return cover


def add_communities(cover, comms):
	cover = copy(cover)
	comm_id = cover.upperBound()
	cover.setUpperBound(comm_id + len(comms))
	for comm in comms:
		for u in comm:
			cover.addToSubset(comm_id, u)
		comm_id += 1
	return cover


# Test how cleaning up the result covers changes the entropy.
def cleanup_test(benchmarks, result_dir):
	for b in benchmarks:
		test_entropy_cleanup(b.get_graph(), b.get_ground_truth(), b.get_cover(),
		                     result_dir)


# Take a cover and try to increase its quality by removing bad nodes/communities.
def test_entropy_cleanup(G, ground_truth, cover, result_dir):
	result = OrderedDict()
	for type in [
		"remove_bad",
		"remove_good",
		"add_good",
		"add_bad"
	]:
		cover_copy = copy(cover)
		communities = [set() for _ in range(cover_copy.upperBound())]
		for u in G.nodes():
			for c in cover_copy.subsetsOf(u):
				communities[c].add(u)
		gt_communities = [set() for _ in range(ground_truth.upperBound())]
		for u in G.nodes():
			for c in ground_truth.subsetsOf(u):
				gt_communities[c].add(u)

		entropy_norm = 100000
		entropy = calc_entropy(G, cover_copy) / entropy_norm
		iterations = 1
		entropy_increased_cnt = 0
		entropy_decreased_cnt = 0
		entropy_unchanged_cnt = 0
		success_cnt = 0
		for _ in range(iterations):
			c_id, community = random_community(communities)
			# c_id = get_worst_community(G, ground_truth, coverCopy)
			# community = communities[c_id]
			_, best_gt_comm = find_best_gt_community(community, gt_communities)

			func = {"remove_bad": remove_bad_node,
			        "remove_good": remove_good_node,
			        "add_good": add_good_node,
			        "add_bad": add_bad_node
			        }[type]
			success = func(G, cover_copy, c_id, community, best_gt_comm)
			if success:
				success_cnt += 1
			# remove_worst_community(G, ground_truth, coverCopy, communities)

			new_entropy = calc_entropy(G, cover_copy) / entropy_norm
			entropy_change = new_entropy - entropy
			if entropy_change > 0.0:
				entropy_increased_cnt += 1
			elif entropy_change < 0.0:
				entropy_decreased_cnt += 1
			elif success:
				entropy_unchanged_cnt += 1
			print(
				"Entropy changed from {old:.5f} to {new:.5f} ({sign}{change:.5f})".format(
					old=entropy, new=new_entropy, sign="+" if entropy_change > 0 else "",
					change=entropy_change))
			entropy = new_entropy

		success_cnt = success_cnt or 1
		result[type] = {"incr": entropy_increased_cnt / success_cnt,
		                "decr": entropy_decreased_cnt / success_cnt,
		                "unchanged": entropy_unchanged_cnt / success_cnt,
		                "success_cnt": success_cnt}
	out_file = open(result_dir + "cleanup.result", 'a')
	output_str = "Random community\n"
	for type, values in result.items():
		s = "+: {incr:.4f}, -: {decr:.4f}, 0: {unchanged:.4f} ({count})\n".format(
			type=type, incr=values["incr"], decr=values["decr"],
			unchanged=values["unchanged"], count=values["success_cnt"])
		output_str += type.ljust(14) + s

	print(output_str)
	out_file.write(output_str)


# Get a list of the nodes of a community, sorted ascending by the number of edges into
# the community.
def worst_nodes_in_comm(G, community):
	nodes = []
	for u in community:
		cnt = internal_edge_cnt(G, community, u)
		nodes.append((cnt / G.degree(u), u))
	nodes.sort()
	return nodes


# Returns the number of edges between the node u and other nodes in the community
def internal_edge_cnt(G, community, u):
	count = 0
	for v in G.neighbors(u):
		if v in community:
			count += 1
	# G.forNeighborsOf(u, lambda v: work(v))
	return count


# Convert the cover (node -> set of its communities) to a list of communities.
# Each community in the list is a list of nodes.
def get_community_vector(G, cover):
	communities = [set() for _ in range(cover.upperBound())]
	for u in G.nodes():
		for c in cover.subsetsOf(u):
			communities[c].add(u)
	return communities


# Calculate the number of edges to other communities.
# Returns a sorted list of (cut_size, community ID) tuples.
def community_cut(G, cover, comm_sizes, community, c_id):
	cuts = defaultdict(lambda: 0)
	for u in community:
		for v in G.neighbors(u):
			for c in cover.subsetsOf(v):
				cuts[c] += 1

	cuts[c_id] = 0
	comm_cuts = sorted([((cut - 1) / (comm_sizes[comm] * comm_sizes[c_id]), comm)  # TODO: What to do with small comms and few edges
	                    for comm, cut in cuts.items()
	                    if cut > 0 and comm_sizes[comm] > 0
	                    ], reverse=True)
	# print(c_id, "Comm size:", comm_sizes[c_id])
	# print([(cut, comm_sizes[comm], cut / (comm_sizes[comm] * comm_sizes[c_id]), comm) for comm, cut in cuts.items() if cut > 0])
	print([(score, cuts[c], comm_sizes[c], c) for score, c in comm_cuts])
	return comm_cuts


# Remove a node from the community that should be in the community
def remove_good_node(G, cover, c_id, community, best_gt_comm):
	good_nodes = best_gt_comm.intersection(community)
	good_node = list(good_nodes)[0]
	cover.removeFromSubset(c_id, good_node)
	community.remove(good_node)
	return True


# Add a node that shouldn't be in the community
def add_bad_node(G, cover, c_id, community, best_gt_comm):
	exclude_nodes = community.union(best_gt_comm)
	while True:
		u = randint(0, G.upperNodeIdBound() - 1)
		if G.hasNode(u) and u not in exclude_nodes:
			break
	cover.addToSubset(c_id, u)
	community.add(u)
	return True


# Remove a node that should not be in the community
def remove_bad_node(G, cover, c_id, community, best_gt_comm):
	bad_nodes = community - best_gt_comm
	if len(bad_nodes) == 0:
		return False
	bad_node = list(bad_nodes)[0]
	cover.removeFromSubset(c_id, bad_node)
	community.remove(bad_node)
	return True


# Add a node that should be in the community
def add_good_node(G, cover, c_id, community, best_gt_comm):
	good_nodes = best_gt_comm - community
	if len(good_nodes) == 0:
		return False
	good_node = list(good_nodes)[0]
	cover.addToSubset(c_id, good_node)
	community.add(good_node)
	return True


def random_community(communities):
	comm = None
	while not comm:
		i = randint(0, len(communities) - 1)
		comm = communities[i]
	assert len(comm) > 0
	return i, comm


# Find the best ground truth community for the given found community.
# Returns the community ID and the node set.
def find_best_gt_community(community, ground_truth_communitites):
	if len(community) == 0:
		return None, set()
	# community = set(community)
	best_score = 0.0
	best_gt_comm_id = none
	for i, gt_comm in enumerate(ground_truth_communitites):
		if len(gt_comm) == 0:
			continue
		overlap = community.intersection(gt_comm)
		if len(overlap) == 0:
			continue

		precision = len(overlap) / len(gt_comm)
		recall = len(overlap) / len(community)
		score = 2 * precision * recall / (precision + recall)

		if score > best_score:
			best_score = score
			best_gt_comm_id = i

	best_comm = set()
	if best_gt_comm_id != none:
		best_comm = ground_truth_communitites[best_gt_comm_id]
	return best_gt_comm_id, best_comm


# Get the community with the lowest F1 score
def get_worst_community(G, ground_truth, cover):
	similarity = CoverF1Similarity(G, cover, ground_truth).run()
	scores = similarity.getValues()
	scores = [(s, i) for i, s in enumerate(scores) if s > 0.0]
	min_score, worst_community = min(scores)
	return worst_community


# Remove the community with the lowest F1 score
def remove_worst_community(G, ground_truth, cover, communities):
	similarity = CoverF1Similarity(G, cover, ground_truth).run()
	scores = similarity.getValues()
	scores = [(s, i) for i, s in enumerate(scores) if s > 0.0]
	min_score, worst_community = min(scores)
	print("Removing community {} with score {}".format(worst_community, min_score))
	for u in communities[worst_community]:
		cover.removeFromSubset(worst_community, u)
	communities[worst_community] = []


# Calculate the F1 score
def calc_f1(graph, cover, refCover):
	similarity = CoverF1Similarity(graph, cover, refCover).run()
	return similarity.getUnweightedAverage()
