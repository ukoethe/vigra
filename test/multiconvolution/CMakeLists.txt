VIGRA_ADD_TEST(test_multiconvolution test.cxx LIBRARIES vigraimpex)

VIGRA_ADD_TEST(test_multiconvolution_speed speedtest.cxx)

VIGRA_ADD_TEST(test_blockwise_convolution blockwise.cxx)

VIGRA_COPY_TEST_DATA(oi_single.gif)
