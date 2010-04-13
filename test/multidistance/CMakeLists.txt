if(HDF5_FOUND)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  
    ADD_DEFINITIONS(${HDF5_CPPFLAGS})

    VIGRA_ADD_TEST(test_multidistance test.cxx LIBRARIES vigraimpex ${HDF5_LIBRARIES})
else()
    VIGRA_ADD_TEST(test_multidistance test.cxx LIBRARIES vigraimpex)
endif()
