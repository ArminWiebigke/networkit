networkit_add_module(community
    AdjustedRandMeasure.cpp
    ClusteringGenerator.cpp
    CommunityDetectionAlgorithm.cpp
    Conductance.cpp
    CoverHubDominance.cpp
    CoverF1Similarity.cpp
    Coverage.cpp
    CutClustering.cpp
    DissimilarityMeasure.cpp
    DynamicNMIDistance.cpp
    EdgeCut.cpp
    GraphClusteringTools.cpp
    GraphStructuralRandMeasure.cpp
    HubDominance.cpp
    IntrapartitionDensity.cpp
    IsolatedInterpartitionConductance.cpp
    IsolatedInterpartitionExpansion.cpp
    JaccardMeasure.cpp
    LPDegreeOrdered.cpp
    LPPotts.cpp
    LocalCoverEvaluation.cpp
	LouvainMapEquation.cpp
    LocalPartitionEvaluation.cpp
    Modularity.cpp
    NMIDistance.cpp
    NodeStructuralRandMeasure.cpp
    OLP.cpp
    PLM.cpp
    PLP.cpp
    ParallelAgglomerativeClusterer.cpp
    PartitionFragmentation.cpp
    PartitionHubDominance.cpp
    PartitionIntersection.cpp
    QualityMeasure.cpp
    SampledGraphStructuralRandMeasure.cpp
    SampledNodeStructuralRandMeasure.cpp
    SLPA.cpp
    StablePartitionNodes.cpp
    cleanup/MergeCommunities.cpp
    cleanup/SignificanceCommunityCleanUp.cpp
    cleanup/SingleCommunityCleanUp.cpp
    cleanup/StochasticDistribution.cpp
	cleanup/SignificanceCalculator.cpp
    egosplitting/EgoNetExtensionAndPartition.cpp
    egosplitting/EgoSplitting.cpp
    egosplitting/ExtendByScore.cpp
    egosplitting/ExtendEgoNetStrategy.cpp
    egosplitting/ExtendSignificance.cpp
    )

networkit_module_link_modules(community
    auxiliary base coarsening components flow graph matching structures)

add_subdirectory(test)
add_subdirectory(cleanup/test)
add_subdirectory(egosplitting/test)

