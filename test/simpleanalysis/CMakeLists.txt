INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})

ADD_EXECUTABLE(test_simpleanalysis test.cxx)

TARGET_LINK_LIBRARIES(test_simpleanalysis vigraimpex ${FFTW3_LIBRARIES})

VIGRA_ADD_TEST(test_simpleanalysis noiseNormalizationTest.xv slantedEdgeMTF.xv lenna128.xv)
