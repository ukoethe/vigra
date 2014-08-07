if(FFTW3_FOUND)
    INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})

    VIGRA_CONFIGURE_THREADING()
      
    VIGRA_ADD_TEST(test_registration test.cxx LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3F_LIBRARIES} ${THREADING_LIBRARIES} vigraimpex)
    VIGRA_COPY_TEST_DATA(nuernberg-1991.png nuernberg-1995.png)
else()
    MESSAGE(STATUS "** WARNING: test_correlation will not be executed")
endif()
