from egosplit.benchmarks.execution.cleanup_functions import merge_overlap_comms, trim_comms, \
	add_communities, remove_small_comms, remove_overlap_comms
from egosplit.external import cleanUpOslom
from egosplit.benchmarks.data_structures.context_timer import ContextTimer


class CleanUpConfig:
	@staticmethod
	def get_clean_up_set(clean_up_set):
		if clean_up_set == 'No Cleanup':
			clean_ups = ['']
		elif clean_up_set == 'all':
			clean_ups = [
				'No Clean Up',
				# 'merge-overl',
				# 'Remove Overlapping',
				'Clean-merge',
				# 'Clean-merge & Remove Overlapping',
				# 'clean-full',
				'Clean-remove',
				'OSLOM-full',
			]
		elif clean_up_set == 'best-ego':
			clean_ups = [
				'No Clean Up',
				'Clean-merge',
			]
		elif clean_up_set == 'best':
			clean_ups = [
				'Clean-new',
			]
		elif clean_up_set == 'new_clean':
			clean_ups = [
				'No Clean Up',
				'Clean-merge',
				'Clean-new',
			]
		elif clean_up_set == 'with-without':
			clean_ups = [
				'No Clean Up',
				'Clean-new',
			]
		elif clean_up_set == 'test':
			clean_ups = [
				'No Clean Up',
				'Clean-merge',
				'Clean-new',
			]
		else:
			raise RuntimeError('No clean-up set provided!')
		if len(clean_ups) > 1:
			clean_ups = ['({:03.0f}){}'.format(i, c) for i, c in enumerate(clean_ups)]
		return clean_ups


class CleanUp:
	def __init__(self, name):
		self.name = name
		self.cover = None
		self.timer = ContextTimer()

	def run(self, graph, cover, ground_truth):
		with self.timer:
			self.cover = clean_up_cover(graph, cover, ground_truth, self.name)

	def get_cover(self):
		return self.cover

	def get_time(self):
		return self.timer.elapsed


def clean_up_cover(graph, cover, ground_truth, clean_up):
	if len(clean_up) > 4 and clean_up[0] == '(' and clean_up[4] == ')':
		clean_up = clean_up[5:]
	if clean_up == '' or clean_up == 'No Clean Up' or clean_up == '*':
		return cover
	if clean_up == 'clean-full':
		cover = clean_up_oslom(graph, cover, bad_groups_strat='remove',
		                       check_unions=True, check_minimality=True,
		                       max_extend=100)
	elif clean_up == 'Clean-new':
		clean_up = SignificanceCommunityCleanUp(graph, cover, 0.1, 0.1, 0.5)
		clean_up.run()
		cover = clean_up.getCover()
	elif clean_up == 'Clean-remove':
		cover = clean_up_oslom(graph, cover)
	elif clean_up == 'clean-shrink':
		cover = clean_up_oslom(graph, cover, cleanup_strategy='check')
	elif clean_up == 'Clean-merge':
		cover = clean_up_oslom(graph, cover, bad_groups_strat='merge')
	elif clean_up == 'Clean-search':
		cover = clean_up_oslom(graph, cover, bad_groups_strat='merge', cleanup_strategy='search')
	elif clean_up == 'Clean-merge & Remove Overlapping':
		cover = remove_overlap_comms(graph, cover, min_overlap=1.0)
		cover = clean_up_oslom(graph, cover, bad_groups_strat='merge')
		cover = remove_overlap_comms(graph, cover, min_overlap=0.7)
	elif clean_up == 'OSLOM-full':
		cover = cleanUpOslom(graph, cover)
	elif clean_up == 'clean-merge-5':
		cover = clean_up_oslom(graph, cover, bad_groups_strat='merge', runs=5)
	elif clean_up == 'clean-merge-shrink':
		cover = clean_up_oslom(graph, cover, bad_groups_strat='merge',
		                       discard_max_extend_groups=False)
	elif clean_up == 'clean-keep':
		cover, bad_groups = clean_up_oslom(graph, cover)
		cover = add_communities(cover, bad_groups)
	elif clean_up == 'trim-gt':
		cover = trim_comms(graph, ground_truth, cover)
	elif clean_up == 'trim-gt,overl':
		cover = trim_comms(graph, ground_truth, cover)
		cover = merge_overlap_comms(graph, cover)
	elif clean_up == 'merge-overl':
		cover = merge_overlap_comms(graph, cover, min_overlap=0.7)
	elif clean_up == 'Remove Overlapping':
		cover = remove_overlap_comms(graph, cover, min_overlap=0.7)
	elif clean_up == 'remove-small':
		pass
	elif clean_up != '':
		raise RuntimeError('\'{}\' is not a valid cleanup option!'.format(clean_up))
	cover = remove_small_comms(graph, cover)
	return cover


def clean_up_oslom(G, cover, threshold=0.1, simple_cleanup=True,
                   runs=1, cleanup_strategy='both', max_extend=2,
                   keep_bad_groups=False, bad_groups_strat='remove',
                   discard_max_extend_groups=True,
                   check_minimality=False, check_unions=False,
                   ):
	bad_groups_filename = 'bad_groups.txt'
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
	if bad_groups_strat == 'remove':
		pass
	elif bad_groups_strat == 'keep':
		args.append('-keep_bad_groups')
	elif bad_groups_strat == 'merge':
		args.append('-merge_discarded')
	if check_minimality:
		args.append('-check_min')
	if discard_max_extend_groups:
		args.append('-discard_max_extend_groups')
	# print(args)
	# args = ['-simple_cleanup',
	#         '-merge_discarded', '-discard_max_extend_groups',
	#         '-max_extend', '2',
	#         '-cup_runs', '1']
	# print(args)
	args = [arg.encode('utf-8') for arg in args]
	# print(args)
	# exit(0)
	# args = []

	cleanAlgo = OslomCleanUp(G, cover, args)
	cleanAlgo.run()
	cleanedCover = cleanAlgo.getCover()

	return cleanedCover
