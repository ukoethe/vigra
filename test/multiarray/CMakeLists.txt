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

VIGRA_ADD_TEST(test_multiarray test.cxx LIBRARIES vigraimpex)

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
    INCLUDE_DIRECTORIES(${SUPPRESS_WARNINGS} ${HDF5_INCLUDE_DIR})
    SET(MULTIARRAY_CHUNKED_LIBRARIES vigraimpex ${HDF5_LIBRARIES} ${MULTIARRAY_CHUNKED_LIBRARIES})
  ELSE(TIFF_FOUND)
    SET(MULTIARRAY_CHUNKED_LIBRARIES vigraimpex ${MULTIARRAY_CHUNKED_LIBRARIES})
  ENDIF()

  VIGRA_ADD_TEST(test_multiarray_chunked test_chunked.cxx
                 LIBRARIES ${MULTIARRAY_CHUNKED_LIBRARIES})
endif()

