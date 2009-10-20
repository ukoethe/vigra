# - Find VIGRANUMPY_DEPENDENCIES
# 
FIND_PACKAGE(PythonInterp)

IF(PYTHONINTERP_FOUND)
    VIGRA_FIND_PACKAGE( Boost COMPONENTS PYTHON )

    FIND_PACKAGE(PythonLibs)

    execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c 
                     "from distutils.sysconfig import *; print get_python_lib()"
                      OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
    FILE(TO_CMAKE_PATH ${PYTHON_SITE_PACKAGES} PYTHON_SITE_PACKAGES)
    SET(VIGRA_NUMPY_INSTALL_DIR ${PYTHON_SITE_PACKAGES}
        CACHE PATH "where to install the vigra Python module")

    execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c 
                     "from numpy.distutils.misc_util import *; print ' '.join(get_numpy_include_dirs())" 
                      RESULT_VARIABLE PYTHON_NUMPY_NOT_FOUND
                      OUTPUT_VARIABLE PYTHON_NUMPY_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

    IF(NOT PYTHON_NUMPY_NOT_FOUND)
        FILE(TO_CMAKE_PATH ${PYTHON_NUMPY_INCLUDE_DIR} PYTHON_NUMPY_INCLUDE_DIR)
    ELSE()
        MESSAGE(STATUS "Could NOT find Python numpy module")
    ENDIF()

    execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c 
                     "import sys; print sys.platform" 
                      OUTPUT_VARIABLE PYTHON_PLATFORM OUTPUT_STRIP_TRAILING_WHITESPACE)

    # handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE if 
    # all listed variables are TRUE
    INCLUDE(FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(VIGRANUMPY_DEPENDENCIES DEFAULT_MSG 
                         PYTHONINTERP_FOUND PYTHONLIBS_FOUND
                         Boost_PYTHON_FOUND PYTHON_NUMPY_INCLUDE_DIR PYTHON_SITE_PACKAGES)

    SET(VIGRANUMPY_INCLUDE_DIR ${PYTHON_INCLUDE_PATH} ${Boost_INCLUDE_DIR} ${PYTHON_NUMPY_INCLUDE_DIR}
        CACHE PATH "include directories needed by VIGRA Python bindings")
    SET(VIGRANUMPY_LIBRARIES ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARIES}
        CACHE FILEPATH "libraries needed by VIGRA Python bindings")
ENDIF()
