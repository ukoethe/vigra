cmake_minimum_required(VERSION 2.6)

##################################################
#
#     basic initializations
#
##################################################

PROJECT(vigra)
SET(vigra_version_short 1.7)
SET(vigra_version 1.7.0)
SET(vigra_year 2010)

SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${vigra_SOURCE_DIR}/config)

include(VigraSetDefaults)

INCLUDE_DIRECTORIES(${vigra_SOURCE_DIR}/include)

IF (MSVC)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_SECURE_SCL=0)
ENDIF ()

##################################################
#
#     search for dependencies
#
##################################################

INCLUDE(VigraFindPackage)
VIGRA_FIND_PACKAGE(TIFF NAMES libtiff_i)
VIGRA_FIND_PACKAGE(JPEG NAMES libjpeg)
VIGRA_FIND_PACKAGE(PNG)
VIGRA_FIND_PACKAGE(FFTW3 NAMES libfftw3-3)
VIGRA_FIND_PACKAGE(HDF5)

FIND_PACKAGE(Doxygen)
FIND_PACKAGE(PythonInterp)

IF(WITH_VIGRANUMPY)
    FIND_PACKAGE( VIGRANUMPY_DEPENDENCIES )
ENDIF()

IF(WITH_VALGRIND)
    FIND_PROGRAM(VALGRIND_EXECUTABLE valgrind)
ENDIF()

##################################################
#
#     setup testing environment
#
##################################################

ADD_CUSTOM_TARGET(check)
ADD_CUSTOM_TARGET(ctest COMMAND ${CMAKE_CTEST_COMMAND})

ADD_CUSTOM_TARGET(check_cpp)
ADD_DEPENDENCIES(check check_cpp)

INCLUDE(CTest)

IF(WITH_VIGRANUMPY)
    ADD_CUSTOM_TARGET(check_python)
    ADD_DEPENDENCIES(check check_python)
ENDIF()

ENABLE_TESTING()

##################################################
#
#     setup documentation generator
#
##################################################

ADD_CUSTOM_TARGET(doc)

##################################################
#
#     setup private experimentation target
#
##################################################

ADD_CUSTOM_TARGET(experiments)

##################################################
#
#     subdirectories
#
##################################################

INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/include)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)
ADD_SUBDIRECTORY(docsrc)

IF(WITH_VIGRANUMPY)
    ADD_SUBDIRECTORY(vigranumpy)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/private)
	ADD_SUBDIRECTORY(private EXCLUDE_FROM_ALL)
ENDIF()


##################################################
#
#     global installation commands
#
##################################################

INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/vigra DESTINATION include)
INSTALL(DIRECTORY ${PROJECT_BINARY_DIR}/include/vigra DESTINATION include)

##################################################
#
#     enable CPack package generation
#
##################################################
INCLUDE(CPackConfig)

# configure version header
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/configVersion.hxx.in
		${PROJECT_BINARY_DIR}/include/vigra/configVersion.hxx
		@ONLY)


##################################################
#
#     enable other CMake projects to find Vigra 
#
##################################################

# config file for CMake FIND_PACKAGE command (for using the build tree)
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/VigraConfigBuild.cmake.in
		${PROJECT_BINARY_DIR}/VigraConfig.cmake
		@ONLY IMMEDIATE)

# config file for CMake FIND_PACKAGE command
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/VigraConfig.cmake.in
		${PROJECT_BINARY_DIR}/lib/vigra/CMake/VigraConfig.cmake
		@ONLY IMMEDIATE)


# export targets:
#export(TARGETS all FILE 
INSTALL(EXPORT vigra-targets DESTINATION lib/vigra)
INSTALL(FILES ${PROJECT_BINARY_DIR}/lib/vigra/CMake/VigraConfig.cmake DESTINATION lib/vigra)
EXPORT(TARGETS vigraimpex FILE vigra-targets.cmake)

##################################################
#
#     status output
#
##################################################

MESSAGE( STATUS "---------------------------------------------------------" )
MESSAGE( STATUS "VIGRA configuration information:" )
MESSAGE( STATUS "---------------------------------------------------------" )

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

IF(WITH_VIGRANUMPY)
    IF(VIGRANUMPY_DEPENDENCIES_FOUND)
        MESSAGE( STATUS "  Using Python libraries: ${VIGRANUMPY_LIBRARIES}" )
    ELSE()
        MESSAGE( STATUS "  Vigranumpy dependencies not found (vigranumpy disabled)" )
    ENDIF()
ENDIF()

MESSAGE( STATUS "---------------------------------------------------------" )
MESSAGE( STATUS "  binaries will be generated in: ${CMAKE_CURRENT_BINARY_DIR}")
IF(DOXYGEN_FOUND AND PYTHONINTERP_FOUND)
	MESSAGE( STATUS "  manuals will be generated in: ${DOCDIR}")
ENDIF()
	MESSAGE( STATUS "---------------------------------------------------------" )

MESSAGE( STATUS "  includes will be installed at: ${CMAKE_INSTALL_PREFIX}/include")
MESSAGE( STATUS "  libraries will be installed at: ${CMAKE_INSTALL_PREFIX}/lib")

IF (MSVC)
  MESSAGE( STATUS "  binaries and DLLs will be installed at: ${CMAKE_INSTALL_PREFIX}/bin")
else()
  MESSAGE( STATUS "  binaries will be installed at: ${CMAKE_INSTALL_PREFIX}/bin")
endif()

IF(DOXYGEN_FOUND AND PYTHONINTERP_FOUND)
  MESSAGE( STATUS "  vigra manuals will be installed at: "
    "${CMAKE_INSTALL_PREFIX}/${DOCINSTALL}/vigra/index.html")
ELSE()
  MESSAGE( STATUS "  vigra manuals cannot be generated (Doxygen or Python not found)")
ENDIF()

IF(WITH_VIGRANUMPY AND VIGRANUMPY_DEPENDENCIES_FOUND)
  MESSAGE( STATUS "  vigranumpy will be installed at ${CMAKE_INSTALL_PREFIX}/${VIGRANUMPY_INSTALL_DIR}" )
  IF(PYTHON_SPHINX)
    MESSAGE( STATUS "  vigranumpy manuals will be installed at: "
      "${CMAKE_INSTALL_PREFIX}/${DOCINSTALL}/vigranumpy/html/index.html")
  ELSE(PYTHON_SPHINX)
    MESSAGE( STATUS "  vigranumpy manuals cannot be generated "
      "(Sphinx missing)")
  ENDIF(PYTHON_SPHINX)
ELSE()
  MESSAGE( STATUS "  vigranumpy will NOT be installed" )
ENDIF()

IF(NOT TIFF_FOUND OR NOT JPEG_FOUND OR NOT PNG_FOUND OR NOT FFTW3_FOUND OR 
   NOT HDF5_FOUND OR NOT VIGRANUMPY_DEPENDENCIES_FOUND)
  MESSAGE( STATUS "" )
  MESSAGE( STATUS "Consider setting DEPENDENCY_SEARCH_PREFIX to find missing libraries" )
ENDIF()

MESSAGE( STATUS "---------------------------------------------------------" )
