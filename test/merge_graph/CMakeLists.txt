if(Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    VIGRA_ADD_TEST(test_adjacency_list_graph_merge test.cxx)
else()
    MESSAGE(WARNING "test_adjacency_list_graph_merge skipped because boost is missing.")
endif()

