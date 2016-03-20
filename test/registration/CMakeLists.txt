VIGRA_CONFIGURE_THREADING()

if(FFTW3_FOUND)
    INCLUDE_DIRECTORIES(${SUPPRESS_WARNINGS} ${FFTW3_INCLUDE_DIR})
    ADD_DEFINITIONS(-DHasFFTW3)

    VIGRA_ADD_TEST(test_registration test.cxx LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3F_LIBRARIES} ${THREADING_LIBRARIES} vigraimpex)
else()
    MESSAGE(STATUS "** WARNING: fftw not found. Fourier-domain registration tests will not be executed")

    VIGRA_ADD_TEST(test_registration test.cxx LIBRARIES ${THREADING_LIBRARIES} vigraimpex)
endif()

VIGRA_COPY_TEST_DATA(nuernberg-1991.png nuernberg-1995.png)
