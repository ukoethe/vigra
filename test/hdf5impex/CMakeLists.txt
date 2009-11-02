if(HDF5_FOUND)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  
    ADD_DEFINITIONS(-DHasHDF5 -D_HDF5USEDLL_ -DHDF5CPP_USEDLL)

    VIGRA_ADD_TEST(test_hdf5impex test.cxx LIBRARIES vigraimpex ${HDF5_LIBRARIES})
else()
    MESSAGE(STATUS "** WARNING: test_hdf5impex will not be executed")
endif()
