
if(Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

    VIGRA_ADD_TEST(test_rag test.cxx)
else()
    MESSAGE(WARNING "test_rag skipped because boost is missing.")
endif()
