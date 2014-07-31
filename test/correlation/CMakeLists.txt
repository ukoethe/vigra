if(FFTW3_FOUND)
    INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})

    VIGRA_ADD_TEST(test_correlation test.cxx LIBRARIES ${FFTW3_LIBRARIES})
else()
    MESSAGE(STATUS "** WARNING: test_correlation will not be executed")
endif()
