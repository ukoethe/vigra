VIGRA_CONFIGURE_THREADING()

if(THREADING_FOUND)
    VIGRA_ADD_TEST(test_threadpool test.cxx LIBRARIES ${THREADING_LIBRARIES})
else()
    MESSAGE(STATUS "** WARNING: No threading implementation found.")
    MESSAGE(STATUS "**          test_threadpool will not be executed on this platform.")
endif()
