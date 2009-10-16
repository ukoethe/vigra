cmake_minimum_required(VERSION 2.6)

PROJECT(vigra)
SET(vigra_version 1.7.0)


ADD_CUSTOM_TARGET(check)

#ENABLE_TESTING()
ADD_CUSTOM_TARGET(test)
ADD_DEPENDENCIES(test check)

if (WIN32 AND ${CMAKE_C_COMPILER} MATCHES gcc.exe)
    SET(IS_CYGWIN 1)
else()
    SET(IS_CYGWIN 0)
endif()

INCLUDE_DIRECTORIES(${vigra_SOURCE_DIR}/include)

SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${vigra_SOURCE_DIR}/config)

SET(CMAKE_INCLUDE_PATH $ENV{HOME}/include
    CACHE PATH "additional include search path (used by Find... macros)")
SET(CMAKE_LIBRARY_PATH $ENV{HOME}/lib
    CACHE PATH "additional library search path (used by Find... macros)")

IF (CMAKE_GENERATOR MATCHES "Visual Studio")
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_SECURE_SCL=0)
    SET(CMAKE_INCLUDE_PATH ${vigra_SOURCE_DIR}/include/external ${CMAKE_INCLUDE_PATH}
        CACHE PATH "additional include search path (used by Find... macros)")
    SET(CMAKE_LIBRARY_PATH ${vigra_SOURCE_DIR}/lib ${CMAKE_LIBRARY_PATH}
        CACHE PATH "additional library search path (used by Find... macros)")
ENDIF (CMAKE_GENERATOR MATCHES "Visual Studio")

SET(TIFF_NAMES libtiff_i)
FIND_PACKAGE(TIFF)

SET(JPEG_NAMES libjpeg)
FIND_PACKAGE(JPEG)

FIND_PACKAGE(PNG)

SET(FFTW3_NAMES libfftw3-3)
FIND_PACKAGE(FFTW3)

FIND_PACKAGE(HDF5)

FIND_PACKAGE(Doxygen)

FIND_PACKAGE(PythonInterp)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)
ADD_SUBDIRECTORY(docsrc)

INSTALL(DIRECTORY ${vigra_SOURCE_DIR}/include/vigra DESTINATION include)
INSTALL(DIRECTORY ${vigra_SOURCE_DIR}/doc/vigra DESTINATION doc)

MESSAGE( STATUS "---------------------------------------------------------" )
MESSAGE( STATUS "VIGRA configuration information:" )
MESSAGE( STATUS "---------------------------------------------------------" )
MESSAGE( STATUS "  includes  will be installed at: ${CMAKE_INSTALL_PREFIX}/include")
MESSAGE( STATUS "  libraries will be installed at: ${CMAKE_INSTALL_PREFIX}/lib")

IF (CMAKE_GENERATOR MATCHES "Visual Studio")
  MESSAGE( STATUS "  binaries and DLLs will be installed at: ${CMAKE_INSTALL_PREFIX}/bin")
else()
  MESSAGE( STATUS "  binaries  will be installed at: ${CMAKE_INSTALL_PREFIX}/bin")
endif()

IF(DOXYGEN_FOUND AND PYTHONINTERP_FOUND)
  MESSAGE( STATUS "  manuals   will be installed at: ${CMAKE_INSTALL_PREFIX}/doc/vigra/index.html")
ELSE()
  MESSAGE( STATUS "  manuals cannot be generated (Doxygen or Python not found)")
ENDIF()

IF(TIFF_FOUND)
  MESSAGE( STATUS "  Using TIFF libraries: ${TIFF_LIBRARIES}" )
ELSE()
  MESSAGE( STATUS "  TIFF libraries not found (TIFF support disabled)" )
ENDIF()

IF(JPEG_FOUND)
  MESSAGE( STATUS "  Using JPEG libraries: ${JPEG_LIBRARIES}" )
ELSE()
  MESSAGE( STATUS "  JPEG libraries not found (JPEG support disabled)" )
ENDIF()

IF(PNG_FOUND)
  MESSAGE( STATUS "  Using PNG  libraries: ${PNG_LIBRARIES}" )
ELSE()
  MESSAGE( STATUS "  PNG libraries not found (PNG support disabled)" )
ENDIF()

IF(TIFF_FOUND)
  MESSAGE( STATUS "  Using FFTW libraries: ${FFTW3_LIBRARIES}" )
ELSE()
  MESSAGE( STATUS "  FFTW libraries not found (FFTW support disabled)" )
ENDIF()

IF(HDF5_FOUND)
  MESSAGE( STATUS "  Using HDF5 libraries: ${HDF5_LIBRARIES}" )
ELSE()
  MESSAGE( STATUS "  HDF5 libraries not found (HDF5 support disabled)" )
ENDIF()

MESSAGE( STATUS "---------------------------------------------------------" )
