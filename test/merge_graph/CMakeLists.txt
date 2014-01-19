if(Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    VIGRA_ADD_TEST(test_rag_merge test.cxx)
else()
    MESSAGE(WARNING "test_rag_merge skipped because boost is missing.")
endif()

