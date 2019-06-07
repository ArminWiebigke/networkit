networkit_add_module(oslom
        LogTable.cpp
        LouvainOslomnet.cpp
        ModuleCollection.cpp
        OslomCleanUp.cpp
        OslomnetEvaluate.cpp
        OslomNetGlobal.cpp
        Parameters.cpp
        StaticNetwork.cpp
        Stochastics.cpp
        WeightedTabdeg.cpp
        utils/Cast.cpp
        utils/Combinatorics.cpp
        utils/DequeNumeric.cpp
        utils/Histograms.cpp
        utils/Print.cpp
        utils/OslomRandom.cpp
        utils/Wsarray.cpp
    )

networkit_module_link_modules(community
    auxiliary base coarsening components flow graph matching structures)

add_subdirectory(test)

