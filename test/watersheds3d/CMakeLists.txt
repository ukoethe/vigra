ADD_EXECUTABLE(test_watersheds3d test.cxx)

TARGET_LINK_LIBRARIES(test_watersheds3d vigraimpex)

VIGRA_ADD_TEST(test_watersheds3d)
