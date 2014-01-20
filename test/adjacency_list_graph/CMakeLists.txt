
if(Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

    VIGRA_ADD_TEST(test_adjacency_list_graph test.cxx)
else()
    MESSAGE(WARNING "test_adjacency_list_graph skipped because boost is missing.")
endif()
