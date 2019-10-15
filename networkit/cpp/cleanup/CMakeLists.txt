networkit_add_module(cleanup
    MergeCommunities.cpp
    SignificanceCommunityCleanUp.cpp
    SingleCommunityCleanUp.cpp
    StochasticDistribution.cpp
    StochasticSignificance.cpp
    )

networkit_module_link_modules(community
    auxiliary base coarsening components flow graph matching structures)

add_subdirectory(test)

