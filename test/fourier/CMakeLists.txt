if(FFTW3_FOUND)
    INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})

    VIGRA_CONFIGURE_THREADING()
      
    VIGRA_ADD_TEST(test_fourier test.cxx LIBRARIES vigraimpex ${FFTW3_LIBRARIES} ${THREADING_LIBRARIES})

    VIGRA_COPY_TEST_DATA(ghouse.gif filter.xv gaborresult.xv)
else()
    MESSAGE(STATUS "** WARNING: test_fourier will not be executed")
endif()
