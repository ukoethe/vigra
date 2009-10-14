INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})

VIGRA_ADD_TEST(test_simpleanalysis test.cxx LIBRARIES vigraimpex ${FFTW3_LIBRARIES})

VIGRA_COPY_TEST_DATA(noiseNormalizationTest.xv slantedEdgeMTF.xv lenna128.xv)
