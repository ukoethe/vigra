ADD_EXECUTABLE(test_volumelabeling test.cxx)

TARGET_LINK_LIBRARIES(test_volumelabeling vigraimpex)

VIGRA_ADD_TEST(test_volumelabeling)
