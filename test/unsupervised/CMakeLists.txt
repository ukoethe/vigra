if(HDF5_FOUND)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  
    ADD_DEFINITIONS(${HDF5_CPPFLAGS})

    VIGRA_ADD_TEST(test_unsupervised test.cxx LIBRARIES vigraimpex ${HDF5_LIBRARIES})
    
    VIGRA_COPY_TEST_DATA(example_data.h5)
else()
    MESSAGE(STATUS "** WARNING: test_unsupervised will not be executed")
endif()