VIGRA_CONFIGURE_THREADING()

if(NOT ${VIGRA_CPP_VERSION})
    message(FATAL_ERROR
            "cmake error: VIGRA_CPP_VERSION not defined yet. "
            "Call VIGRA_DETECT_CPP_VERSION() from the main CMakeLists file.")
endif()

# c++11 'auto' is used in the test as well as in multi_chunked.hxx
string(COMPARE LESS ${VIGRA_CPP_VERSION} "201103" NO_CXX11)
if(NO_CXX11)
    MESSAGE(STATUS "** WARNING: You are compiling in c++98 mode.")
    MESSAGE(STATUS "**          Chunked convolution tests will be skipped.")
    MESSAGE(STATUS "**          Add -std=c++11 to CMAKE_CXX_FLAGS to enable this tests.")
else()
    if(NOT THREADING_FOUND)
        MESSAGE(STATUS "** WARNING: No threading implementation found.")
        MESSAGE(STATUS "**          Chunked convolution tests will be skipped.")
    else()
        VIGRA_ADD_TEST(test_chunked_convolution test_chunked_convolution.cxx
            LIBRARIES vigraimpex ${THREADING_LIBRARIES})
    endif()
endif()
