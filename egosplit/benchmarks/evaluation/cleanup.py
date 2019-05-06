from random import randint
from collections import OrderedDict

from networkit.community import CoverF1Similarity
from networkit.structures import Cover
from networkit import none
from egosplit.external import calc_entropy


# Clean up the result covers.
def cleanup(benchmarks, result_dir):
	for b in benchmarks:
		test_entropy_cleanup(b.get_graph(), b.get_ground_truth(), b.get_cover(), result_dir)


# Take a cover and try to increase its quality by removing bad nodes/communities.
def test_entropy_cleanup(G, ground_truth, cover, result_dir):
	result = OrderedDict()
	for type in [
		"remove_bad",
		"remove_good",
		"add_good",
		"add_bad"
	]:
		coverCopy = cover  # TODO: copy
		communities = [set() for _ in range(coverCopy.upperBound())]
		for u in G.nodes():
			for c in coverCopy.subsetsOf(u):
				communities[c].add(u)
		gt_communities = [set() for _ in range(ground_truth.upperBound())]
		for u in G.nodes():
			for c in ground_truth.subsetsOf(u):
				gt_communities[c].add(u)

		entropy_norm = 100000
		entropy = calc_entropy(G, coverCopy) / entropy_norm
		iterations = 100
		entropy_increased_cnt = 0
		entropy_decreased_cnt = 0
		entropy_unchanged_cnt = 0
		success_cnt = 0
		for _ in range(iterations):
			c_id, community = random_community(communities)
			# c_id = get_worst_community(G, ground_truth, coverCopy)
			# community = communities[c_id]
			best_gt_comm = find_best_gt_community(community, gt_communities)

			func = {"remove_bad": remove_bad_node,
			        "remove_good": remove_good_node,
			        "add_good": add_good_node,
			        "add_bad": add_bad_node
			        }[type]
			success = func(G, coverCopy, c_id, community, best_gt_comm)
			if success:
				success_cnt += 1
			# remove_worst_community(G, ground_truth, coverCopy, communities)

			new_entropy = calc_entropy(G, coverCopy) / entropy_norm
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


# Find the best ground truth community for the given found community
def find_best_gt_community(community, ground_truth_communitites):
	if len(community) == 0:
		return none
	# community = set(community)
	best_score = 0.0
	best_gt_comm = none
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
			best_gt_comm = i

	return ground_truth_communitites[best_gt_comm]


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
