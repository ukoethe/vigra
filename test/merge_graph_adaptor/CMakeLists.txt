if(Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    VIGRA_ADD_TEST(test_merge_graph_adaptor test.cxx)
else()
    MESSAGE(WARNING "test_merge_graph_adaptor skipped because boost is missing.")
endif()

