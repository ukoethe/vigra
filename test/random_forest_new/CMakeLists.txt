if(HDF5_FOUND)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
    ADD_DEFINITIONS(${HDF5_CPPFLAGS} -DHasHDF5)
    VIGRA_ADD_TEST(test_random_forest_new test.cxx LIBRARIES vigraimpex ${HDF5_LIBRARIES})
else()
    VIGRA_ADD_TEST(test_random_forest_new test.cxx)
endif()

add_subdirectory(data)

