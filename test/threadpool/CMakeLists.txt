VIGRA_CONFIGURE_THREADING()

if(THREADING_FOUND)
    VIGRA_ADD_TEST(test_threadpool test.cxx)
else()
    MESSAGE(STATUS "** WARNING: Your compiler does not support std::thread.")
    MESSAGE(STATUS "**          test_threadpool will not be executed on this platform.")
endif()
