if(FFTW3_FOUND)
    INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})
    ADD_DEFINITIONS(-DHasFFTW3)
    
    VIGRA_ADD_TEST(test_simpleanalysis test.cxx LIBRARIES vigraimpex ${FFTW3_LIBRARIES})
else()
    VIGRA_ADD_TEST(test_simpleanalysis test.cxx LIBRARIES vigraimpex)
endif()

VIGRA_COPY_TEST_DATA(noiseNormalizationTest.xv slantedEdgeMTF.xv lenna128.xv)
