if(FFTW3_FOUND)
    INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})

    VIGRA_ADD_TEST(test_features test.cxx LIBRARIES vigraimpex ${FFTW3_LIBRARIES})

    # VIGRA_COPY_TEST_DATA(ghouse.gif)
else()
    MESSAGE(STATUS "** WARNING: test_features will not be executed")
endif()

