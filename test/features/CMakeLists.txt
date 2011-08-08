if(FFTW3_FOUND)
    INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})
    
    IF(HDF5_FOUND)
        ADD_DEFINITIONS(-DHasHDF5 ${HDF5_CPPFLAGS})
        INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
    ENDIF(HDF5_FOUND)

    VIGRA_ADD_TEST(test_features test.cxx LIBRARIES vigraimpex ${FFTW3_LIBRARIES} ${FFTW3F_LIBRARIES})

    # VIGRA_COPY_TEST_DATA(ghouse.gif)
else()
    MESSAGE(STATUS "** WARNING: test_features will not be executed")
endif()

