from collections import OrderedDict

from egosplit.benchmarks.data_structures.algorithms import GroundTruth, GceAlgorithm, MosesAlgorithm, OslomAlgorithm, \
	DemonAlgorithm, EgoSplitAlgorithm
from egosplit.benchmarks.execution.cleanup import CleanUpConfig
from egosplit.external import partitionLeiden, partitionInfomap
from networkit.community import PLM, PLP, LPPotts, LouvainMapEquationFactory, PLMFactory, PLPFactory
# from egosplit.benchmarks.data_structures.algorithms import *


class EgoSplitClusteringAlgorithmsConfig:
	@staticmethod
	def get(ego_part_algos):
		partition_algos = OrderedDict()

		if ego_part_algos == 'local' or ego_part_algos == 'global':
			partition_algos['PLP'] = [PLPFactory(1, 20)]
			partition_algos['PLM'] = [PLMFactory(True, 1.0, 'none')]
			partition_algos['Potts'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition()]
			partition_algos['Infomap'] = [lambda g: partitionInfomap(g)]
			partition_algos['Surprise'] = [lambda g: partitionLeiden(g, 'surprise')]
			partition_algos['Leiden'] = [lambda g: partitionLeiden(g, 'modularity')]
			partition_algos['LM-Map'] = [LouvainMapEquationFactory(True)]

			new_p_algos = {}
			for name, p_algos in partition_algos.items():
				if ego_part_algos == 'local':
					new_p_algos[name + ' + LM-Map'] = [p_algos[0],
					                                   LouvainMapEquationFactory(True)]
				if ego_part_algos == 'global':
					new_p_algos['Leiden + ' + name] = [lambda g: partitionLeiden(g, 'modularity'),
					                                   p_algos[0]]
			partition_algos = new_p_algos

		if ego_part_algos == 'best':
			partition_algos['Leiden + LM-Map'] = [lambda g: partitionLeiden(g, 'modularity'),
			                                      LouvainMapEquationFactory(True)]
		if ego_part_algos == 'two_best':
			partition_algos['Leiden + LM-Map'] = [lambda g: partitionLeiden(g, 'modularity'),
			                                      LouvainMapEquationFactory(True)]
			partition_algos['LM-Map + LM-Map'] = [lambda g: partitionInfomap(g),
			                                         lambda g: partitionLeiden(g, 'surprise')]
		if ego_part_algos == 'Leiden/Infomap + Infomap':
			partition_algos['Leiden + Infomap'] = [lambda g: partitionLeiden(g, 'modularity'),
			                                       lambda g: partitionInfomap(g)]
			partition_algos['Infomap + Infomap'] = [lambda g: partitionInfomap(g),
			                                        lambda g: partitionInfomap(g)]
		if ego_part_algos == 'leiden local':
			partition_algos['Leiden + Infomap'] = [lambda g: partitionLeiden(g, 'modularity'),
			                                       lambda g: partitionInfomap(g)]
		if ego_part_algos == 'fast':
			partition_algos['PLP + PLM'] = [lambda g: PLP(g, 1, 20).run().getPartition(),
			                                lambda g: PLM(g, True, 1.0,
			                                              'none').run().getPartition()]
		if ego_part_algos == 'test':
			partition_algos['PLP + LM-Map'] = [PLPFactory(1, 20), LouvainMapEquationFactory(True)]
			# def error_func(g):
			# 	raise RuntimeError("Expected error")
			# partition_algos['error'] = [error_func, error_func]
			# partition_algos['default'] = []
			partition_algos['PLM (fac)'] = [PLMFactory(True, 1.0, 'none')]
			partition_algos['PLM'] = [lambda g: PLM(g, True, 1.0, 'none').run().getPartition()]
			partition_algos['PLM + PLM (both fac)'] = [PLMFactory(True, 1.0, 'none'),
			                                           PLMFactory(True, 1.0, 'none')]
		# partition_algos['PLM + PLM'] = [lambda g: PLM(g, True, 1.0, 'none').run().getPartition(),
		#                                 lambda g: PLM(g, True, 1.0, 'none').run().getPartition()]
		# partition_algos['PLM(fac) + PLM'] = [PLMFactory(True, 1.0, 'none'),
		#                                      lambda g: PLM(g, True, 1.0, 'none').run().getPartition()]
		# partition_algos['PLM + PLM(fac)'] = [lambda g: PLM(g, True, 1.0, 'none').run().getPartition(),
		#                                      PLMFactory(True, 1.0, 'none')]

		return partition_algos


class EgoSplitParameterConfig:
	@staticmethod
	def get(ego_parameter_config):
		ego_parameters = OrderedDict()
		standard = {
			'Extend EgoNet Strategy': 'None',
			'Extend and Partition Iterations': 1,
			'partitionFromGroundTruth': 'No',
			'numEgoNetsStored': 2000,
			'Cleanup': 'No',
			'connectPersonas': 'Yes',
			'normalizePersonaCut': 'No',
			'connectPersonasStrat': 'spanning',
			'maxPersonaEdges': 1,
			'normalizePersonaWeights': 'unweighted',
			'iterationWeight': 'No',
		}
		extend_standard = {
			**standard,
			'Extend and Partition Iterations': 1,
			'Maximum Extend Factor': 5,
			'addNodesExponent': 0.5,
			'minNodeDegree': 2,
			'keepOnlyTriangles': 'No',
		}
		edge_scores_standard = {
			**extend_standard,
			'Extend EgoNet Strategy': 'Edges',
			'Edges Score Strategy': 'Edges pow 2 div Degree',
		}
		significance_scores_standard = {
			**extend_standard,
			'Extend EgoNet Strategy': 'Significance',
			'Significance Base Extend': 'None',
			'maxSignificance': 0.1,
			'sortGroups': 'Significance',
			'orderedStatPos': 0.0,
			'useSigMemo': 'No',
			'minEdgesToGroupSig': 1,
			'maxGroupsConsider': 99,
			'secondarySigExtRounds': 99,
			'signMerge': 'Yes',
			'Extend and Partition Iterations': 3,
			'onlyCheckSignOfMaxCandidates': 'Yes',
			'Check Candidates Factor': 10,
			'onlyUpdatedCandidates': 'Yes',
		}

		if 'no-extend' in ego_parameter_config:
			ego_parameters['No Extension'] = standard
		if 'edges-score' in ego_parameter_config:
			for score in ['Edges', 'Edges div Degree', 'Edges pow 2 div Degree', 'Random', 'Significance']:
				name = 'Extend: {}'.format(score)
				ego_parameters[name] = {
					**edge_scores_standard,
					'Edges Score Strategy': score,
				}
		if 'edges-significance' in ego_parameter_config:
			for score in ['Edges pow 2 div Degree', 'Random', 'Significance']:
				name = 'Extend: {}'.format(score)
				ego_parameters[name] = {
					**edge_scores_standard,
					'Edges Score Strategy': score,
				}
		if 'edges-factor' in ego_parameter_config:
			for factor in [1, 2, 3, 5, 10, 20]:
				name = r'$\alpha = {}$'.format(factor)
				ego_parameters[name] = {
					**edge_scores_standard,
					'Maximum Extend Factor': factor,
				}
		if 'sig-merge' in ego_parameter_config:
			for merge in [False, True]:
				name = 'Single + Merged Clusters' if merge else 'Single Clusters'
				ego_parameters[name] = {
					**significance_scores_standard,
					'signMerge': 'Yes' if merge else 'No',
					'onlyCheckSignOfMaxCandidates': 'No',
					'secondarySigExtRounds': 0,
					'Extend and Partition Iterations': 1,
				}
		if 'sig-max-candidates' in ego_parameter_config:
			for max_factor in [1, 2, 3, 5, 10, 20, 10000]:
				name = r'$\gamma = {}$'.format(max_factor)
				if max_factor == 10000:
					name = 'All candidates'
				ego_parameters[name] = {
					**significance_scores_standard,
					'onlyCheckSignOfMaxCandidates': 'Yes',
					'Check Candidates Factor': max_factor,
					'secondarySigExtRounds': 0,
					'Extend and Partition Iterations': 1,
				}
		if 'sig-ext-iter' in ego_parameter_config:
			for iterations in [0, 1, 2, 3, 5, 10, 100]:
				name = '{} Iteration{}'.format(iterations, '' if iterations == 1 else 's')
				ego_parameters[name] = {
					**significance_scores_standard,
					# 'onlyCheckSignOfMaxCandidates': 'No',
					'secondarySigExtRounds': iterations,
					'onlyUpdatedCandidates': 'No',
					'Extend and Partition Iterations': 1,
				}
		if 'sig-check-updated' in ego_parameter_config:
			for updated in [False, True]:
				name = 'Only Improved' if updated else 'All'
				ego_parameters[name] = {
					**significance_scores_standard,
					'onlyUpdatedCandidates': 'Yes' if updated else 'No',
					# 'onlyCheckSignOfMaxCandidates': 'No',
					'Extend and Partition Iterations': 1,
				}
		if 'sig-cluster-iter' in ego_parameter_config:
			for iterations in [1, 2, 3, 5, 8]:
				name = '$I_c$ = {}'.format(iterations)
				ego_parameters[name] = {
					**significance_scores_standard,
					'Extend and Partition Iterations': iterations,
				}
		if 'sig-mem' in ego_parameter_config:
			for memoize in [True, False]:
				name = 'Memoize' if memoize else 'Calculate'
				ego_parameters[name] = {
					**significance_scores_standard,
					'useSigMemo': 'Yes' if memoize else 'No',
				}
		if 'extend' in ego_parameter_config:
			ego_parameters['EdgesScore'] = {
				**edge_scores_standard,
			}
			ego_parameters['Significance'] = {
				**significance_scores_standard,
			}
		if 'connect-persona' in ego_parameter_config:
			ego_parameters['EdgesScore | No Connection'] = {
				**edge_scores_standard,
				'connectPersonas': 'No',
			}
			ego_parameters['EdgesScore | Max Spanning Unweighted'] = {
				**edge_scores_standard,
				'connectPersonas': 'Yes',
				'connectPersonasStrat': 'spanning',
				'normalizePersonaCut': 'No',
				'maxPersonaEdges': 1,
				'normalizePersonaWeights': 'unweighted',
				'iterationWeight': 'No',
			}
			ego_parameters['EdgesScore | All Density Max Weight 1'] = {
				**edge_scores_standard,
				'connectPersonas': 'Yes',
				'connectPersonasStrat': 'all',
				'normalizePersonaCut': 'density',
				'normalizePersonaWeights': 'max1',
			}
			ego_parameters['EdgesScore | All Unweighted'] = {
				**edge_scores_standard,
				'connectPersonas': 'Yes',
				'connectPersonasStrat': 'all',
				'normalizePersonaCut': 'No',
				'normalizePersonaWeights': 'unweighted',
			}
		if 'edges' in ego_parameter_config:
			ego_parameters['EdgesScore'] = {
				**edge_scores_standard,
			}
		if 'steps' in ego_parameter_config:
			ego_parameters['Base'] = {
				**standard,
				'connectPersonas': 'No',
			}
			ego_parameters['EdgesScore'] = {
				**edge_scores_standard,
				'connectPersonas': 'No',
			}
			ego_parameters['E + Spanning'] = {
				**edge_scores_standard,
			}
			ego_parameters['E + S + Cleanup'] = {
				**edge_scores_standard,
				'Cleanup': 'Yes',
			}
		if 'best' in ego_parameter_config:
			ego_parameters[''] = {
				**edge_scores_standard,
				'Cleanup': 'Yes',
			}

		if 'test' in ego_parameter_config:
			ego_parameters['EdgesScore'] = {
				**edge_scores_standard,
			}

		return ego_parameters


class OtherAlgorithms:
	@staticmethod
	def get(algo_list):
		algos = []  # Elements are tuples (AlgorithmObject, list of clean up procedures)
		algos.append((GroundTruth(), ['']))
		if not algo_list:
			return algos

		algorithm_create = {
			'GCE': lambda: GceAlgorithm(alpha=1.1),
			'Moses': lambda: MosesAlgorithm(),
			'Oslom': lambda: OslomAlgorithm(),
			'Demon': lambda: DemonAlgorithm(),
			'Ego-original': lambda: EgoSplitAlgorithm(
				'Ego-original', original_ego_parameters(),
				lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition(),
				lambda g: LPPotts(g, 0, 1, 20, True).run().getPartition()
			),
			'Ego-base': lambda: EgoSplitAlgorithm(
				'Ego-base', original_ego_parameters(),
				lambda g: partitionLeiden(g, 'modularity'),
				lambda g: partitionInfomap(g)
			),
			'Ego-optimized': lambda: EgoSplitAlgorithm(
				'Ego-optimized', optimized_ego_parameters(),
				lambda g: partitionLeiden(g, 'modularity'),
				LouvainMapEquationFactory(True)
			),
		}
		for algo in algo_list:
			algos.append((algorithm_create[algo](), ['']))

		# algos.append(OlpAlgorithm())
		# algos.append((GceAlgorithm('GCE', alpha=1.0), ['', 'OSLOM-merge']))
		# algos.append(GceAlgorithm('GCE_1.0_clean', alpha=1.0, clean_up='OSLOM-merge'))
		# algos.append((MosesAlgorithm(), ['', 'OSLOM-merge']))

		return algos



def get_ego_algos(ego_part_algos, ego_parameter_config, clean_up_set, store_ego_nets):
	if not ego_part_algos or not ego_parameter_config:
		return []
	algos = []
	part_algos = EgoSplitClusteringAlgorithmsConfig.get(ego_part_algos)

	ego_parameters = EgoSplitParameterConfig.get(ego_parameter_config)
	for parameter_set in ego_parameters.values():
		parameter_set['storeEgoNet'] = 'Yes' if store_ego_nets else 'No'
	clean_ups = CleanUpConfig.get_clean_up_set(clean_up_set)
	algos += create_egosplit_algorithms(part_algos, ego_parameters, clean_ups)

	return algos


def egosplit_partition_algorithms(ego_part_algos):
	partition_algos = OrderedDict()

	if ego_part_algos == 'local' or ego_part_algos == 'global':
		partition_algos['PLP'] = [lambda g: PLP(g, 1, 20).run().getPartition()]
		partition_algos['PLM'] = [lambda g: PLM(g, True, 1.0, 'none').run().getPartition()]
		# partition_algos['Potts_0.01'] = [lambda g: LPPotts(g, 0.01, 1, 20).run().getPartition()]
		# partition_algos['Potts_0.05'] = [lambda g: LPPotts(g, 0.05, 1, 20).run().getPartition()]
		partition_algos['Potts'] = [lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition()]
		# partition_algos['LPPotts_par'] = [
		# 	lambda g: LPPotts(g, 0.1, 1, 20).run().getPartition(),
		# 	lambda g: LPPotts(g, 0, 1, 20, True).run().getPartition()]
		partition_algos['Infomap'] = [lambda g: partitionInfomap(g)]
		partition_algos['Surprise'] = [lambda g: partitionLeiden(g, 'surprise')]
		partition_algos['Leiden'] = [lambda g: partitionLeiden(g, 'modularity')]
		# partition_algos['Leiden_Sig'] = [lambda g: partitionLeiden(g, 'significance')]
		# partition_algos['Leiden_SigRes'] = [lambda g: leidenSignificance(g)]

		new_p_algos = {}
		for name, p_algos in partition_algos.items():
			if ego_part_algos == 'local':
				new_p_algos[name + ' + Infomap'] = [p_algos[0],
				                                    lambda g: partitionInfomap(g)]  # Infomap global
			if ego_part_algos == 'global':
				new_p_algos['Infomap + ' + name] = [lambda g: partitionInfomap(g), p_algos[0]]
				new_p_algos['Leiden + ' + name] = [lambda g: partitionLeiden(g, 'modularity'),
				                                   p_algos[0]]
		partition_algos = new_p_algos

	if ego_part_algos == 'best':
		partition_algos['Leiden + Infomap'] = [lambda g: partitionLeiden(g, 'modularity'),
		                                       lambda g: partitionInfomap(g)]
		partition_algos['Infomap + Surprise'] = [lambda g: partitionInfomap(g),
		                                         lambda g: partitionLeiden(g, 'surprise')]
	if ego_part_algos == 'Leiden/Infomap + Infomap':
		partition_algos['Leiden + Infomap'] = [lambda g: partitionLeiden(g, 'modularity'),
		                                       lambda g: partitionInfomap(g)]
		partition_algos['Infomap + Infomap'] = [lambda g: partitionInfomap(g),
		                                        lambda g: partitionInfomap(g)]
	if ego_part_algos == 'standard':
		partition_algos['Leiden + Infomap'] = [lambda g: partitionLeiden(g, 'modularity'),
		                                       lambda g: partitionInfomap(g)]
	if ego_part_algos == 'leiden local':
		partition_algos['Leiden + Infomap'] = [lambda g: partitionLeiden(g, 'modularity'),
		                                       lambda g: partitionInfomap(g)]
	if ego_part_algos == 'fast':
		partition_algos['PLP + PLM'] = [lambda g: PLP(g, 1, 20).run().getPartition(),
		                                lambda g: PLM(g, True, 1.0, 'none').run().getPartition()]

	return partition_algos


def original_ego_parameters():
	return {
		'weightFactor': 0,
		'weightOffset': 1,
		'storeEgoNet': 'No',
		'Extend EgoNet Strategy': 'None',
		'Extend and Partition Iterations': 1,
		'Maximum Extend Factor': 0,
		'addNodesExponent': 0,
		'partitionFromGroundTruth': 'No',
		'connectPersonas': 'No',
		'numEgoNetsStored': 2000,
	}


def optimized_ego_parameters():
	return {
		'partitionFromGroundTruth': 'No',
		'numEgoNetsStored': 2000,
		'connectPersonas': 'Yes',
		'normalizePersonaCut': 'No',
		'connectPersonasStrat': 'spanning',
		'maxPersonaEdges': 1,
		'normalizePersonaWeights': 'unweighted',
		'iterationWeight': 'No',
		#
		'Extend EgoNet Strategy': 'Edges',
		'Edges Score Strategy': 'Edges pow 2 div Degree',
		'Extend and Partition Iterations': 1,
		'Maximum Extend Factor': 5,
		'addNodesExponent': 0.5,
		'minNodeDegree': 2,
		'keepOnlyTriangles': 'No',
		#
		'Cleanup': 'Yes',
	}


def create_egosplit_algorithms(partition_algos, ego_parameters, clean_ups):
	algos = []
	i = 0
	for part_name in partition_algos:
		for para_name, parameters in ego_parameters.items():
			name = '{}{}{}'.format('Ego({:03.0f})'.format(i), part_name,
			                       ' | ' + para_name if para_name else '')
			algo = EgoSplitAlgorithm(
				name,
				parameters,
				*partition_algos[part_name])
			algos.append((algo, clean_ups))
			i += 1
	return algos
