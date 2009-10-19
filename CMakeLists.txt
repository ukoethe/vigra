cmake_minimum_required(VERSION 2.6)

PROJECT(vigra)
SET(vigra_version 1.7.0)

ADD_CUSTOM_TARGET(check)

#ENABLE_TESTING()
ADD_CUSTOM_TARGET(test)
ADD_DEPENDENCIES(test check)

if (WIN32 AND ${CMAKE_C_COMPILER} MATCHES /bin/gcc.exe)
    SET(IS_CYGWIN 1)
else()
    SET(IS_CYGWIN 0)
endif()

INCLUDE_DIRECTORIES(${vigra_SOURCE_DIR}/include)

SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${vigra_SOURCE_DIR}/config)

SET(DEPENDENCY_SEARCH_PREFIX ""
    CACHE PATH "additional search prefix (used by Find... macros)")

IF (CMAKE_GENERATOR MATCHES "Visual Studio")
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_SECURE_SCL=0)
ENDIF (CMAKE_GENERATOR MATCHES "Visual Studio")

INCLUDE(VigraFindPackage)
VIGRA_FIND_PACKAGE(TIFF libtiff_i)
VIGRA_FIND_PACKAGE(JPEG libjpeg)
VIGRA_FIND_PACKAGE(PNG)
VIGRA_FIND_PACKAGE(FFTW3 libfftw3-3)
VIGRA_FIND_PACKAGE(HDF5)

FIND_PACKAGE(Doxygen)
FIND_PACKAGE(PythonInterp)

IF(PYTHONINTERP_FOUND)
    SET(BOOST_INCLUDEDIR $ENV{HOME}/include/boost_1_34_1)
    SET(BOOST_LIBRARYDIR $ENV{HOME}/lib)
    FIND_PACKAGE( Boost COMPONENTS python )
    
    MESSAGE( STATUS "Boost_FOUND: ${Boost_FOUND}" )
    MESSAGE( STATUS "Using BOOST_PYTHON libraries: ${Boost_PYTHON_LIBRARIES}" )
    FIND_PACKAGE(PythonLibs)
    
    execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c 
                     "from distutils.sysconfig import *; print get_python_lib()"
                      OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
    FILE(TO_CMAKE_PATH ${PYTHON_SITE_PACKAGES} PYTHON_SITE_PACKAGES)
    SET(VIGRA_NUMPY_INSTALL_DIR ${PYTHON_SITE_PACKAGES}/vigra)

    execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c 
                     "from numpy.distutils.misc_util import *; print ' '.join(get_numpy_include_dirs())" 
                      OUTPUT_VARIABLE PYTHON_NUMPY_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    FILE(TO_CMAKE_PATH ${PYTHON_NUMPY_INCLUDE_DIR} PYTHON_NUMPY_INCLUDE_DIR)

    execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c 
                     "import sys; print sys.platform" 
                      OUTPUT_VARIABLE PYTHON_PLATFORM OUTPUT_STRIP_TRAILING_WHITESPACE)
ENDIF()

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(python)
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

IF(FFTW3_FOUND)
  MESSAGE( STATUS "  Using FFTW libraries: ${FFTW3_LIBRARIES}" )
ELSE()
  MESSAGE( STATUS "  FFTW libraries not found (FFTW support disabled)" )
ENDIF()

IF(HDF5_FOUND)
  MESSAGE( STATUS "  Using HDF5 libraries: ${HDF5_LIBRARIES}" )
ELSE()
  MESSAGE( STATUS "  HDF5 libraries not found (HDF5 support disabled)" )
ENDIF()

IF(NOT TIFF_FOUND OR NOT JPEG_FOUND OR NOT PNG_FOUND OR NOT FFTW3_FOUND OR NOT HDF5_FOUND)
  MESSAGE( STATUS "" )
  MESSAGE( STATUS "Consider setting DEPENDENCY_SEARCH_PREFIX to find missing libraries" )
ENDIF()

MESSAGE( STATUS "---------------------------------------------------------" )
