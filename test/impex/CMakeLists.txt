ADD_EXECUTABLE(test_impex test.cxx)

TARGET_LINK_LIBRARIES(test_impex vigraimpex)

VIGRA_ADD_TEST(test_impex lenna.xv lenna_gifref.xv lennafloat.xv lennafloatrgb.xv lennargb.xv no-image.txt)
