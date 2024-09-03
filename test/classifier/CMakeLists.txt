if(HDF5_FOUND)
    INCLUDE_DIRECTORIES(${SUPPRESS_WARNINGS} ${HDF5_INCLUDE_DIR})
    ADD_DEFINITIONS(${HDF5_CPPFLAGS} -DHasHDF5)

    VIGRA_ADD_TEST(test_classifier test.cxx LIBRARIES vigraimpex ${HDF5_LIBRARIES})
    VIGRA_ADD_TEST(classifier_speed_comparison speed_comparison.cxx LIBRARIES ${HDF5_LIBRARIES})
else()
    MESSAGE(STATUS "** WARNING: test_classifier::RFHDF5Test() will not be executed")

    VIGRA_ADD_TEST(test_classifier test.cxx )
    VIGRA_ADD_TEST(classifier_speed_comparison speed_comparison.cxx)
endif()

add_subdirectory(data)

