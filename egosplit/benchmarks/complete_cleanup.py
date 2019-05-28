from .cleanup import merge_overlap_comms, merge_comms_entropy, trim_comms, add_communities, \
	remove_small_comms
from ..external import cleanUpOSLOM


def clean_up_cover(graph, cover, ground_truth, clean_up):
	if clean_up == "":
		return cover
	if clean_up == "OSLOM_full":
		cover, _ = cleanUpOSLOM(graph, cover, clean_bad="merge", check_minimality=True,
		                        max_extend=100)
	elif clean_up == "OSLOM_remove":
		cover, _ = cleanUpOSLOM(graph, cover)
	elif clean_up == "OSLOM_shrink":
		cover, _ = cleanUpOSLOM(graph, cover, cleanup_strategy="check")
	elif clean_up == "OSLOM_merge":
		cover, _ = cleanUpOSLOM(graph, cover, clean_bad="merge")
	elif clean_up == "OSLOM_keep":
		cover, bad_groups = cleanUpOSLOM(graph, cover)
		cover = add_communities(cover, bad_groups)
	elif clean_up == "OSLOM_keep-merge_E":
		cover, bad_groups = cleanUpOSLOM(graph, cover)
		cover = merge_overlap_comms(graph, cover)
		cover = merge_comms_entropy(graph, cover, bad_groups)
		cover = merge_overlap_comms(graph, cover)
		# cover = remove_comms_entropy(graph, cover)
	elif clean_up == "trim_gt":
		cover = trim_comms(graph, ground_truth, cover)
	elif clean_up == "trim_gt-overl":
		cover = trim_comms(graph, ground_truth, cover)
		cover = merge_overlap_comms(graph, cover)
	elif clean_up == "trim_gt-merge_E":
		cover = trim_comms(graph, ground_truth, cover)
		cover = merge_overlap_comms(graph, cover)
		cover = merge_comms_entropy(graph, cover)
		cover = merge_overlap_comms(graph, cover)
	elif clean_up != "":
		raise RuntimeError("\"{}\" is not a valid cleanup option!".format(clean_up))
	cover = remove_small_comms(graph, cover)
	return cover