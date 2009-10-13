ADD_EXECUTABLE(test_tensorimaging test.cxx)

TARGET_LINK_LIBRARIES(test_tensorimaging vigraimpex)

VIGRA_ADD_TEST(test_tensorimaging l2.xv riesz00.xv riesz10.xv riesz01.xv riesz20.xv riesz11.xv riesz02.xv boundaryTensor.xv l2_boundary1.xv l2_boundary.xv l2_hourglass.xv l2_get.xv)
