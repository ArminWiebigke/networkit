networkit_add_test(community CommunityGTest
    auxiliary generators io overlap structures)

networkit_add_test(community MapEquationGTest
	auxiliary generators io overlap structures)

networkit_add_benchmark(community CommunityDetectionBenchmark
	auxiliary centrality io structures)

networkit_add_benchmark(community MapEquationBenchmark
	auxiliary centrality io structures)
