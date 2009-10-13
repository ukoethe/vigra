ADD_EXECUTABLE(test_multidistance test.cxx)

TARGET_LINK_LIBRARIES(test_multidistance vigraimpex)

VIGRA_ADD_TEST(test_multidistance)
