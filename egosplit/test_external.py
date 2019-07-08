from external import *
from networkit.stopwatch import clockit


@clockit
def test():
	G, C = genLFR()
	for u in G.nodes():
		# if u > 10:
		# 	continue
		egoGraph = G.subgraphFromNodes(G.neighbors(u))
		# partitionLeiden(egoGraph, 'surprise')
		partition = leidenSignificance(egoGraph)
		# print(partition)

test()