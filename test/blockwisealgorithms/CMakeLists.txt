VIGRA_CONFIGURE_THREADING()

if(THREADING_FOUND)
    # VIGRA_ADD_TEST(test_blockwiselabeling test_labeling.cxx LIBRARIES ${THREADING_LIBRARIES}) # FIXME
    VIGRA_ADD_TEST(test_blockwisewatersheds test_watersheds.cxx LIBRARIES ${THREADING_LIBRARIES})
    VIGRA_ADD_TEST(test_blockwiseconvolution test_convolution.cxx LIBRARIES ${THREADING_LIBRARIES})
else()
    MESSAGE(STATUS "** WARNING: No threading implementation found.")
    MESSAGE(STATUS "**          test_blockwiselabeling will not be executed on this platform.")
    MESSAGE(STATUS "**          test_blockwisewatersheds will not be executed on this platform.")
    MESSAGE(STATUS "**          test_blockwiseconvolution will not be executed on this platform.")
endif()