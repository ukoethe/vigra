VIGRA_ADD_TEST(test_multiconvolution test.cxx LIBRARIES vigraimpex)

VIGRA_ADD_TEST(test_multiconvolution_speed speedtest.cxx)

VIGRA_COPY_TEST_DATA(oi_single.gif lenna.gif deriche_2_0_10.xv
  deriche_2_1_10.xv  deriche_2_2_10.xv  deriche_3_0_10.xv
  deriche_3_1_10.xv  deriche_3_2_10.xv  deriche_4_0_10.xv
  deriche_4_1_10.xv  deriche_4_2_10.xv
  deriche_2_0_15.xv  deriche_2_1_15.xv  deriche_2_2_15.xv
  deriche_3_0_15.xv  deriche_3_1_15.xv  deriche_3_2_15.xv
  deriche_4_0_15.xv  deriche_4_1_15.xv  deriche_4_2_15.xv
  deriche_2_0_5.xv   deriche_2_1_5.xv   deriche_2_2_5.xv
  deriche_3_0_5.xv   deriche_3_1_5.xv   deriche_3_2_5.xv deriche_4_0_5.xv   deriche_4_1_5.xv   deriche_4_2_5.xv)
