INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})

ADD_EXECUTABLE(test_fourier test.cxx)

TARGET_LINK_LIBRARIES(test_fourier vigraimpex ${FFTW3_LIBRARIES})

VIGRA_ADD_TEST(test_fourier ghouse.gif filter.xv gaborresult.xv)
