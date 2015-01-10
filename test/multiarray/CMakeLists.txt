IF(JPEG_FOUND)
  ADD_DEFINITIONS(-DHasJPEG)
ENDIF(JPEG_FOUND)

IF(PNG_FOUND)
  ADD_DEFINITIONS(-DHasPNG)
ENDIF(PNG_FOUND)

IF(TIFF_FOUND)
  ADD_DEFINITIONS(-DHasTIFF)
ENDIF(TIFF_FOUND)

FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/impex)

# Check cpp version
if(NOT ${VIGRA_CPP_VERSION})
    message(FATAL_ERROR 
            "cmake error: VIGRA_CPP_VERSION not defined yet. "
            "Call VIGRA_DETECT_CPP_VERSION() from the main CMakeLists file." )
endif()

# multiarray/test.cxx uses 'auto' from c++11.
string(COMPARE LESS ${VIGRA_CPP_VERSION} "201103" NO_CXX11)
if(NO_CXX11)
    MESSAGE(STATUS "** WARNING: You are compiling in C++98 mode.")
    MESSAGE(STATUS "**          Multiarray tests will be skipped.")
    MESSAGE(STATUS "**          Add -std=c++11 to CMAKE_CXX_FLAGS to enable multiarray tests.")
else()
    VIGRA_ADD_TEST(test_multiarray test.cxx LIBRARIES vigraimpex)
endif()    

# Even with C++11, a working threading implementation is needed for running multiarray_chunked tests.
VIGRA_CONFIGURE_THREADING()

if(NOT THREADING_FOUND)
    MESSAGE(STATUS "** WARNING: Your compiler does not support C++ threading.")
    MESSAGE(STATUS "**          test_multiarray_chunked will not be executed on this platform.")
    if(NOT WITH_BOOST_THREAD)
        MESSAGE(STATUS "**          Try to run cmake with '-DWITH_BOOST_THREAD=1' to use boost threading.")
    endif()
else()
  SET(MULTIARRAY_CHUNKED_LIBRARIES ${THREADING_LIBRARIES})

  IF(HDF5_FOUND)
    ADD_DEFINITIONS(-DHasHDF5 ${HDF5_CPPFLAGS})
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
    SET(MULTIARRAY_CHUNKED_LIBRARIES vigraimpex ${HDF5_LIBRARIES} ${MULTIARRAY_CHUNKED_LIBRARIES})
  ELSE(TIFF_FOUND)
    SET(MULTIARRAY_CHUNKED_LIBRARIES vigraimpex ${MULTIARRAY_CHUNKED_LIBRARIES})
  ENDIF()

  VIGRA_ADD_TEST(test_multiarray_chunked test_chunked.cxx 
                 LIBRARIES ${MULTIARRAY_CHUNKED_LIBRARIES})
endif()

