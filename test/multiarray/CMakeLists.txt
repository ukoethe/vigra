ADD_EXECUTABLE(test_multiarray test.cxx)

TARGET_LINK_LIBRARIES(test_multiarray vigraimpex)

VIGRA_ADD_TEST(test_multiarray)
