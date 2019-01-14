from egosplit.external import genLFR
from networkit.community import EgoSplitting, CoverF1Similarity, PLM, PLP


def simpleTest(graph, refCover, algo, localClustering=None, globalClustering=None):
    algo = EgoSplitting(graph, localClustering, globalClustering)
    cover = algo.run().getCover()
    similarity = CoverF1Similarity(graph, cover, refCover)
    similarity.run()
    return similarity


def startBenchmarks():
    graphNames = ['0.01', '0.1', '0.3']
    graphArgs = {
        '0.01': {'N': 1000, 'k': 25, 'maxk': 50, 'mu': 0.01, 'minc': 20, 'maxc': 50, 'on': 500, 'om': 3},
        '0.1': {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.1, 'minc': 5, 'maxc': 50, 'on': 100, 'om': 2},
        '0.3': {'N': 1000, 'k': 10, 'maxk': 50, 'mu': 0.3, 'minc': 5, 'maxc': 50, 'on': 100, 'om': 2},
    }
    partitionAlgos = {
        'PLP': lambda g: PLP(g).run().getPartition(),
        'PLM': lambda g: PLM(g).run().getPartition()
    }
    iterations = 10

    graphs = {}
    for graphName in graphNames:
        args = graphArgs[graphName]
        graphs[graphName] = []
        for i in range(0, iterations):
            graphs[graphName].append(genLFR(**args))

    for graphName in graphNames:
        print("Graph", graphName)
        for algoName in partitionAlgos.keys():
            algo = partitionAlgos[algoName]
            print("\t" + algoName)
            for graph, cover in graphs[graphName]:
                similarity = simpleTest(graph, cover, algo)
                print("\t\t", similarity.getUnweightedAverage())


startBenchmarks()