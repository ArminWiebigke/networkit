from collections import OrderedDict

from egosplit.benchmarks.data_structures.algorithms import GroundTruth, GceAlgorithm, MosesAlgorithm, OslomAlgorithm, \
	DemonAlgorithm, EgoSplitAlgorithm
from egosplit.benchmarks.execution.cleanup import CleanUpConfig
from egosplit.external import partitionLeiden, partitionInfomap
from networkit.community import PLM, PLP, LPPotts, LouvainMapEquationFactory, PLMFactory, PLPFactory, LPPottsFactory


class OtherAlgorithms:
	@staticmethod
	def get(algo_list):
		algos = []  # Elements are tuples (AlgorithmObject, list of clean up procedures)
		algos.append((GroundTruth(), ['']))
		if not algo_list:
			return algos

		algorithm_construction = {
			'GCE': lambda: GceAlgorithm(alpha=1.1),
			'Moses': lambda: MosesAlgorithm(),
			'Oslom': lambda: OslomAlgorithm(),
			'Demon': lambda: DemonAlgorithm(),
			'Ego-original': lambda: EgoSplitAlgorithm(
				'Ego-original',
				original_ego_parameters(),
				LPPottsFactory(0.1, 1, 20, False),
				LPPottsFactory(0, 1, 20, True)
			),
			'Ego-optimized': lambda: EgoSplitAlgorithm(
				'Ego-optimized', {}
			),
		}
		for algo in algo_list:
			algos.append((algorithm_construction[algo](), ['']))

		return algos


def original_ego_parameters():
	return {
		'Extend EgoNet Strategy': 'None',
		'connectPersonas': 'No',
		'Cleanup': 'No',
	}

