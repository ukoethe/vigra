ADD_EXECUTABLE(test_image test.cxx)

TARGET_LINK_LIBRARIES(test_image vigraimpex)

VIGRA_ADD_TEST(test_image lenna.xv)
