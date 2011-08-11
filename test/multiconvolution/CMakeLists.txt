VIGRA_ADD_TEST(test_multiconvolution test.cxx)

VIGRA_ADD_TEST(test_multiconvolution_speed speedtest.cxx)

VIGRA_ADD_TEST(test_multiconvolution_scaled test_scaled.cxx LIBRARIES vigraimpex)

VIGRA_COPY_TEST_DATA(oi_single.png)
