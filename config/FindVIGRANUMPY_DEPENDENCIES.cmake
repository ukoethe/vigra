# - Find VIGRANUMPY_DEPENDENCIES
# 
FIND_PACKAGE(PythonInterp)

IF(PYTHONINTERP_FOUND)
    VIGRA_FIND_PACKAGE( Boost COMPONENTS python )

    FIND_PACKAGE(PythonLibs)
    IF(NOT PYTHONLIBS_FOUND)
        # fallback when standard search does not work
        execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c 
                         "import sys; skip = 2 if sys.platform.startswith('win') else 1; print 'python' + sys.version[0:3:skip]"
                          OUTPUT_VARIABLE PYTHON_LIBRARY_NAME OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c 
                         "import sys; print sys.exec_prefix"
                          OUTPUT_VARIABLE PYTHON_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
        FIND_LIBRARY(PYTHON_LIBRARIES ${PYTHON_LIBRARY_NAME} "${PYTHON_PREFIX}/libs")
        execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c
                        "from distutils.sysconfig import *; print get_python_inc()"
                         OUTPUT_VARIABLE PYTHON_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
        SET(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE}
            CACHE PATH "Path to Python include files"
            FORCE)
        IF(PYTHON_LIBRARIES AND PYTHON_INCLUDE_PATH)
           SET(PYTHONLIBS_FOUND TRUE)
        ENDIF()
    ENDIF()

    ######################################################################
    #
    #      find default install directory for Python modules
    #      (usually PYTHONDIR/Lib/site-packages)
    #
    ######################################################################
    IF(NOT DEFINED VIGRANUMPY_INSTALL_DIR OR VIGRANUMPY_INSTALL_DIR MATCHES "^$")
        execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c 
                         "from distutils.sysconfig import *; print get_python_lib()"
                          OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
        FILE(TO_CMAKE_PATH ${PYTHON_SITE_PACKAGES} VIGRANUMPY_INSTALL_DIR)
    ENDIF()
    SET(VIGRANUMPY_INSTALL_DIR ${VIGRANUMPY_INSTALL_DIR}
        CACHE PATH "where to install the VIGRA Python package" FORCE)

    ######################################################################
    #
    #      find numpy include directory
    #      (usually below PYTHONDIR/Lib/site-packages/numpy)
    #
    ######################################################################
    execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c 
                     "from numpy.distutils.misc_util import *; print ' '.join(get_numpy_include_dirs())" 
                      RESULT_VARIABLE PYTHON_NUMPY_NOT_FOUND
                      OUTPUT_VARIABLE PYTHON_NUMPY_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

    IF(NOT PYTHON_NUMPY_NOT_FOUND)
        MESSAGE(STATUS "Searching for Python numpy: ok")
        FILE(TO_CMAKE_PATH ${PYTHON_NUMPY_INCLUDE_DIR} PYTHON_NUMPY_INCLUDE_DIR)
    ELSE()
        MESSAGE(STATUS "Could NOT find Python numpy ('import numpy.distutils.misc_util' failed)")
    ENDIF()

    ######################################################################
    #
    #      check if nosetests are installed
    #
    ######################################################################
    execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c 
                     "import nose" 
                      RESULT_VARIABLE PYTHON_NOSETESTS_NOT_FOUND)

    IF(NOT PYTHON_NOSETESTS_NOT_FOUND)
        MESSAGE(STATUS "Searching for Python nosetests: ok")
    ELSE()
        MESSAGE(STATUS "Could NOT find Python nosetests ('import nose' failed)")
    ENDIF()

    ######################################################################
    #
    #      check if sphinx documentation generator is installed
    #
    ######################################################################
    find_program ( PYTHON_SPHINX sphinx-build "${PYTHON_PREFIX}/Scripts")

    IF(NOT PYTHON_SPHINX)
        MESSAGE(STATUS "Could NOT find sphinx documentation generator")
    ELSE()
        MESSAGE(STATUS "Searching for sphinx documentation generator: ok")
    ENDIF()

    ######################################################################
    #
    #      find Python platform
    #
    ######################################################################
    execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c 
                     "import sys; print sys.platform" 
                      OUTPUT_VARIABLE PYTHON_PLATFORM OUTPUT_STRIP_TRAILING_WHITESPACE)

    ######################################################################
    #
    #      set outputs
    #
    ######################################################################
    INCLUDE(FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(VIGRANUMPY_DEPENDENCIES DEFAULT_MSG 
                         PYTHONINTERP_FOUND PYTHONLIBS_FOUND
                         Boost_PYTHON_FOUND PYTHON_NUMPY_INCLUDE_DIR VIGRANUMPY_INSTALL_DIR)

    IF(NOT VIGRANUMPY_INCLUDE_DIRS OR VIGRANUMPY_INCLUDE_DIRS MATCHES "-NOTFOUND")
        SET(VIGRANUMPY_INCLUDE_DIRS ${PYTHON_INCLUDE_PATH} ${Boost_INCLUDE_DIR} ${PYTHON_NUMPY_INCLUDE_DIR})
    ENDIF()    
    SET(VIGRANUMPY_INCLUDE_DIRS ${VIGRANUMPY_INCLUDE_DIRS}
        CACHE PATH "include directories needed by VIGRA Python bindings"
        FORCE)
    IF(NOT VIGRANUMPY_LIBRARIES OR VIGRANUMPY_LIBRARIES MATCHES "-NOTFOUND")
        SET(VIGRANUMPY_LIBRARIES ${PYTHON_LIBRARIES} ${Boost_PYTHON_LIBRARIES})
    ENDIF()    
    SET(VIGRANUMPY_LIBRARIES ${VIGRANUMPY_LIBRARIES}
        CACHE FILEPATH "libraries needed by VIGRA Python bindings"
        FORCE)
ENDIF()
