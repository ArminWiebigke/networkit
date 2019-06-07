from copy import copy

from .cleanup import merge_overlap_comms, merge_comms_entropy, trim_comms, \
	add_communities, \
	remove_small_comms
from networkit.community import OslomCleanUp


def cleanUpOslom(G, cover, threshold=0.1, simple_cleanup=True,
                 runs=1, cleanup_strategy='both', max_extend=2,
                 keep_bad_groups=False, bad_groups_strat="remove",
                 discard_max_extend_groups=True,
                 check_minimality=False, check_unions=False,
                 ):
	bad_groups_filename = "bad_groups.txt"
	args = [
		'-cup_runs', str(runs),
		'-cu_strat', cleanup_strategy,
		'-max_extend', str(max_extend),
		'-threshold', str(threshold),
		'-bad_groups_file', bad_groups_filename,
	]

	if check_unions:
		args.append('-check_unions')
	if keep_bad_groups:
		args.append('-keep_bad')
	if simple_cleanup:
		args.append('-simple_cleanup')
	if bad_groups_strat == "remove":
		pass
	elif bad_groups_strat == "keep":
		args.append('-keep_bad_groups')
	elif bad_groups_strat == "merge":
		args.append('-merge_discarded')
	if check_minimality:
		args.append('-check_min')
	if discard_max_extend_groups:
		args.append('-discard_max_extend_groups')

	args = [arg.encode('utf-8') for arg in args]

	cleanAlgo = OslomCleanUp(copy(G), copy(cover), args)
	cleanAlgo.run()
	cleanedCover = cleanAlgo.getCover()

	bad_groups_file = open(bad_groups_filename, 'r')
	bad_groups = []
	for group in bad_groups_file:
		nodes = group.split(' ')[:-1]  # Last item is the newline character
		bad_groups.append([int(u) for u in nodes])
	return cleanedCover, bad_groups


def clean_up_cover(graph, cover, ground_truth, clean_up):
	if clean_up == "":
		return cover
	if clean_up == "Oslom-full":
		cover, _ = cleanUpOslom(graph, cover, bad_groups_strat="remove",
		                        check_unions=True, check_minimality=True,
		                        max_extend=100)
	elif clean_up == "Oslom-remove":
		cover, _ = cleanUpOslom(graph, cover)
	elif clean_up == "Oslom-shrink":
		cover, _ = cleanUpOslom(graph, cover, cleanup_strategy="check")
	elif clean_up == "Oslom-merge":
		cover, _ = cleanUpOslom(graph, cover, bad_groups_strat="merge")
	elif clean_up == "Oslom-merge-shrink":
		cover, _ = cleanUpOslom(graph, cover, bad_groups_strat="merge",
		                        discard_max_extend_groups=False)
	elif clean_up == "Oslom-keep":
		cover, bad_groups = cleanUpOslom(graph, cover)
		cover = add_communities(cover, bad_groups)
	elif clean_up == "Oslom-keep,merge-E":
		cover, bad_groups = cleanUpOslom(graph, cover)
		cover = merge_overlap_comms(graph, cover)
		cover = merge_comms_entropy(graph, cover, bad_groups)
		cover = merge_overlap_comms(graph, cover)
	# cover = remove_comms_entropy(graph, cover)
	elif clean_up == "trim-gt":
		cover = trim_comms(graph, ground_truth, cover)
	elif clean_up == "trim-gt,overl":
		cover = trim_comms(graph, ground_truth, cover)
		cover = merge_overlap_comms(graph, cover)
	elif clean_up == "trim-gt,merge-E":
		cover = trim_comms(graph, ground_truth, cover)
		cover = merge_overlap_comms(graph, cover)
		cover = merge_comms_entropy(graph, cover)
		cover = merge_overlap_comms(graph, cover)
	elif clean_up == "remove-small":
		pass
	elif clean_up != "":
		raise RuntimeError("\"{}\" is not a valid cleanup option!".format(clean_up))
	cover = remove_small_comms(graph, cover)
	return cover
