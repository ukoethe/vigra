cmake_minimum_required(VERSION 2.6)

##################################################
#
#     basic initializations
#
##################################################

PROJECT(vigra)

# read the current version from configVersion.hxx
file(READ ${vigra_SOURCE_DIR}/include/vigra/configVersion.hxx VIGRA_VERSION_FILE)
string(REGEX MATCH "VIGRA_VERSION_MAJOR[ \t\n]+[^ \t\n]+" VIGRA_VERSION_MAJOR ${VIGRA_VERSION_FILE})
string(REGEX REPLACE "VIGRA_VERSION_MAJOR[ \t\n]" "" VIGRA_VERSION_MAJOR ${VIGRA_VERSION_MAJOR})
string(REGEX MATCH "VIGRA_VERSION_MINOR[ \t\n]+[^ \t\n]+" VIGRA_VERSION_MINOR ${VIGRA_VERSION_FILE})
string(REGEX REPLACE "VIGRA_VERSION_MINOR[ \t\n]" "" VIGRA_VERSION_MINOR ${VIGRA_VERSION_MINOR})
string(REGEX MATCH "VIGRA_VERSION_PATCH[ \t\n]+[^ \t\n]+" VIGRA_VERSION_PATCH ${VIGRA_VERSION_FILE})
string(REGEX REPLACE "VIGRA_VERSION_PATCH[ \t\n]" "" VIGRA_VERSION_PATCH ${VIGRA_VERSION_PATCH})

SET(vigra_version_short ${VIGRA_VERSION_MAJOR}.${VIGRA_VERSION_MINOR})
SET(vigra_version ${VIGRA_VERSION_MAJOR}.${VIGRA_VERSION_MINOR}.${VIGRA_VERSION_PATCH})
message(STATUS "Configuring VIGRA version ${vigra_version}")

SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${vigra_SOURCE_DIR}/config)

include(VigraSetDefaults)
include(VigraCMakeUtils)

INCLUDE_DIRECTORIES(${vigra_SOURCE_DIR}/include)

IF(VIGRA_STATIC_LIB)
    SET(LIBTYPE STATIC)
    ADD_DEFINITIONS(-DVIGRA_STATIC_LIB)
ELSE(VIGRA_STATIC_LIB)
    SET(LIBTYPE SHARED)
ENDIF(VIGRA_STATIC_LIB)

IF (MSVC)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_DEPRECATE)
ENDIF ()

# if(CMAKE_COMPILER_IS_GNUCXX)
  # exec_program(
      # ${CMAKE_CXX_COMPILER}
      # ARGS                    --version
      # OUTPUT_VARIABLE _compiler_output)
  # string(REGEX REPLACE ".*([0-9]\\.[0-9]\\.[0-9]).*" "\\1"  gcc_compiler_version ${_compiler_output})
  # if(gcc_compiler_version VERSION_LESS "4.4.0")
      # message(FATAL_ERROR "VIGRA requires at least gcc 4.4")
  # else()
      # message(STATUS "Using gcc ${gcc_compiler_version} [${CMAKE_CXX_COMPILER}]")
  # endif()
# endif()

##################################################
#
#     search for dependencies
#
##################################################

INCLUDE(VigraFindPackage)
VIGRA_FIND_PACKAGE(TIFF NAMES libtiff)
VIGRA_FIND_PACKAGE(JPEG NAMES libjpeg)
VIGRA_FIND_PACKAGE(PNG)
VIGRA_FIND_PACKAGE(FFTW3 NAMES libfftw3-3 libfftw-3.3)
VIGRA_FIND_PACKAGE(FFTW3F NAMES libfftw3f-3 libfftwf-3.3)
IF(WITH_OPENEXR)
    VIGRA_FIND_PACKAGE(OpenEXR)
ENDIF()
IF(WITH_HDF5)
    VIGRA_FIND_PACKAGE(HDF5)
ENDIF()
IF(WITH_BOOST_GRAPH)
    IF(WITH_VIGRANUMPY)
        VIGRA_FIND_PACKAGE( Boost 1.40.0 COMPONENTS python )
    else()
        VIGRA_FIND_PACKAGE( Boost 1.40.0 )
    endif()
ENDIF()
IF(WITH_LEMON)
    VIGRA_FIND_PACKAGE(LEMON)
ENDIF()

SET(DOXYGEN_SKIP_DOT TRUE)
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


OPTION(CREATE_CTEST_TARGETS "Create special targets for CTest" OFF)
IF(CREATE_CTEST_TARGETS)
  INCLUDE(CTest)
ENDIF()

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

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)
ADD_SUBDIRECTORY(docsrc)

IF(WITH_VIGRANUMPY)
    ADD_SUBDIRECTORY(vigranumpy)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/private/CMakeLists.txt)
    ADD_SUBDIRECTORY(private EXCLUDE_FROM_ALL)
ENDIF()


set(LIB_SUFFIX "" CACHE STRING "Define suffix of lib directory name (32/64)" )

##################################################
#
#     global installation commands
#
##################################################

INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/vigra
        DESTINATION include)

##################################################
#
#     enable CPack package generation
#
##################################################
INCLUDE(CPackConfig)

##################################################
#
#     add target vigra-$version-src.tar.gz
#
##################################################

find_program(GIT_EXECUTABLE 
              NAMES git git.exe git.cmd
              HINTS $ENV{ProgramFiles}/Git/bin
              DOC "git command line client")

CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/package-src.cmake.in
                ${PROJECT_BINARY_DIR}/package-src.cmake
                @ONLY IMMEDIATE)

add_custom_target(PACKAGE_SRC_TAR
                   COMMAND ${CMAKE_COMMAND} -P package-src.cmake
                   WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
                   COMMENT "Creating ${PROJECT_BINARY_DIR}/vigra-${vigra_version}-src.tar.gz")        

ADD_DEPENDENCIES(PACKAGE_SRC_TAR check)
ADD_DEPENDENCIES(PACKAGE_SRC_TAR doc_cpp)
IF(WITH_VIGRANUMPY AND PYTHON_SPHINX)
    ADD_DEPENDENCIES(PACKAGE_SRC_TAR doc_python)
ENDIF()

##################################################
#
#  enable other non-CMake projects to find Vigra
#
##################################################

IF(NOT WIN32)
# vigra-config script for (non-CMake-based) third party packages
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/vigra-config.in
        ${PROJECT_BINARY_DIR}/bin/vigra-config
        @ONLY IMMEDIATE)
INSTALL(FILES ${PROJECT_BINARY_DIR}/bin/vigra-config
        DESTINATION bin
        PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE)
ENDIF(NOT WIN32)

##################################################
#
#     enable other CMake projects to find Vigra
#
##################################################

# config file for CMake FIND_PACKAGE command (for using the build tree)
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/VigraConfigBuild.cmake.in
        ${PROJECT_BINARY_DIR}/VigraConfig.cmake
        @ONLY IMMEDIATE)
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/VigraConfigVersion.cmake.in
        ${PROJECT_BINARY_DIR}/VigraConfigVersion.cmake
        @ONLY IMMEDIATE)

# config file for CMake FIND_PACKAGE command
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/VigraConfig.cmake.in
        ${PROJECT_BINARY_DIR}/lib/vigra/CMake/VigraConfig.cmake
        @ONLY IMMEDIATE)
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/VigraConfigVersion.cmake.in
        ${PROJECT_BINARY_DIR}/lib/vigra/CMake/VigraConfigVersion.cmake
        @ONLY IMMEDIATE)

# export targets:
INSTALL(EXPORT vigra-targets
        DESTINATION lib${LIB_SUFFIX}/vigra)
INSTALL(FILES ${PROJECT_BINARY_DIR}/lib/vigra/CMake/VigraConfig.cmake
            ${PROJECT_BINARY_DIR}/lib/vigra/CMake/VigraConfigVersion.cmake
        DESTINATION lib${LIB_SUFFIX}/vigra)
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

IF(OPENEXR_FOUND)
    MESSAGE( STATUS "  Using OpenEXR  libraries: ${OPENEXR_LIBRARIES}" )
ELSEIF(NOT WITH_OPENEXR)
    MESSAGE( STATUS "  OpenEXR disabled by user (WITH_OPENEXR=0)" )
ELSE()
    MESSAGE( STATUS "  OpenEXR libraries not found (OpenEXR support disabled)" )
ENDIF()

IF(FFTW3_FOUND)
    MESSAGE( STATUS "  Using FFTW libraries: ${FFTW3_LIBRARIES}" )
ELSE()
    MESSAGE( STATUS "  FFTW libraries not found (FFTW support disabled)" )
ENDIF()

IF(HDF5_FOUND)
    MESSAGE( STATUS "  Using HDF5 libraries: ${HDF5_LIBRARIES}" )
ELSEIF(NOT WITH_HDF5)
    MESSAGE( STATUS "  HDF5 disabled by user (WITH_HDF5=0)" )
ELSE()
    MESSAGE( STATUS "  HDF5 libraries not found (HDF5 support disabled)" )
ENDIF()

IF(Boost_FOUND)
    MESSAGE( STATUS "  Using Boost Graph Library: ${Boost_INCLUDE_DIR}/boost/graph" )
ELSEIF(NOT WITH_BOOST_GRAPH)
    MESSAGE( STATUS "  Boost Graph Library disabled by user (WITH_BOOST_GRAPH=0)" )
ELSE()
    MESSAGE( STATUS "  Boost Graph Library not found (support disabled)" )
ENDIF()

IF(LEMON_FOUND)
    MESSAGE( STATUS "  Using LEMON graph library: ${LEMON_LIBRARY}" )
ELSEIF(NOT WITH_LEMON)
    MESSAGE( STATUS "  LEMON graph library disabled by user (WITH_LEMON=0)" )
ELSE()
    MESSAGE( STATUS "  LEMON graph library not found (support disabled)" )
ENDIF()

IF(WITH_VIGRANUMPY)
    IF(VIGRANUMPY_DEPENDENCIES_FOUND)
        MESSAGE( STATUS "  Using Python libraries: ${VIGRANUMPY_LIBRARIES}" )
        MESSAGE( STATUS "  Using Numpy includes:   ${PYTHON_NUMPY_INCLUDE_DIR}" )
    ELSE()
        MESSAGE( STATUS "  Vigranumpy dependencies not found (vigranumpy disabled)" )
    ENDIF()
ELSE()
  MESSAGE( STATUS "  Vigranumpy disabled by user (WITH_VIGRANUMPY=0)" )
ENDIF()

MESSAGE( STATUS "---------------------------------------------------------" )

IF(VIGRA_STATIC_LIB)
    MESSAGE( STATUS "  building static lib" )
ELSE()
    MESSAGE( STATUS "  building shared lib" )
ENDIF()

MESSAGE( STATUS "  binaries will be generated in: ${CMAKE_CURRENT_BINARY_DIR}")
IF(DOXYGEN_FOUND AND PYTHONINTERP_FOUND)
    MESSAGE( STATUS "  manuals will be generated in: ${DOCDIR}")
ENDIF()

MESSAGE( STATUS "---------------------------------------------------------" )

MESSAGE( STATUS "  includes will be installed at: ${CMAKE_INSTALL_PREFIX}/include")
MESSAGE( STATUS "  libraries will be installed at: ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")

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
