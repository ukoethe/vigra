

if(Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

    VIGRA_ADD_TEST(test_graph_algorithm test.cxx)
else()
    MESSAGE(WARNING "test_graph_algorithm skipped because boost is missing.")
endif()
