# - Find VIGRANUMPY_DEPENDENCIES
#
MESSAGE(STATUS "Checking VIGRANUMPY_DEPENDENCIES")

IF(NOT PYTHONINTERP_FOUND)
    FIND_PACKAGE(PythonInterp ${PYTHON_VERSION})
ENDIF()

IF(PYTHONINTERP_FOUND)

    # Note:
    #  'FIND_PACKAGE(PythonLibs)' is unreliable because results are often inconsistent
    #  with the Python interpreter found previously (e.g. libraries or includes
    #  from incompatible installations). Thus, we ask Python itself for the information.
    #

    ######################################################################
    #
    #      find Python prefix
    #
    ######################################################################
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                     "import sys; print(sys.exec_prefix)"
                      OUTPUT_VARIABLE PYTHON_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                     "import sys; print(sys.version.split()[0])"
                      OUTPUT_VARIABLE PYTHON_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                     "import sys; print(sys.version_info[0])"
                      OUTPUT_VARIABLE PYTHON_VERSION_MAJOR OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                     "import sys; print(sys.version_info[1])"
                      OUTPUT_VARIABLE PYTHON_VERSION_MINOR OUTPUT_STRIP_TRAILING_WHITESPACE)

    MESSAGE(STATUS "Using Python ${PYTHON_VERSION} at ${PYTHON_EXECUTABLE}")

    ######################################################################
    #
    #      find Python includes
    #
    ######################################################################
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                    "from distutils.sysconfig import *; print(get_python_inc())"
                    OUTPUT_VARIABLE PYTHON_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
    SET(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE}
        CACHE PATH "Path to Python include files"
        FORCE)

    IF(PYTHON_INCLUDE_PATH)
        MESSAGE(STATUS "Found Python includes:  ${PYTHON_INCLUDE_PATH}")
    ELSE()
        MESSAGE(STATUS "Could NOT find Python includes")
    ENDIF()

    ######################################################################
    #
    #      find Python library
    #
    ######################################################################
    IF(APPLE AND ${PYTHON_PREFIX} MATCHES ".*framework.*")
        SET(PYTHON_LIBRARIES "${PYTHON_PREFIX}/Python"
            CACHE FILEPATH "Python libraries"
            FORCE)
    ELSE()
        IF(WIN32)
            set(PYTHON_LIBRARY_NAME python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
        ELSE()
            execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                             "from distutils.sysconfig import *; print(get_config_var('LDLIBRARY'))"
                              OUTPUT_VARIABLE PYTHON_LIBRARY_NAME OUTPUT_STRIP_TRAILING_WHITESPACE)
            execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                             "from distutils.sysconfig import *; print(get_config_var('LIBDIR'))"
        			           OUTPUT_VARIABLE PYTHON_LIBRARY_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
        ENDIF()
        FIND_LIBRARY(PYTHON_LIBRARIES ${PYTHON_LIBRARY_NAME} HINTS "${PYTHON_LIBRARY_PREFIX}" "${PYTHON_PREFIX}"
                     PATH_SUFFIXES lib lib64 libs DOC "Python libraries")
        unset(PYTHON_LIBRARY_PREFIX)
    ENDIF()

    IF(PYTHON_LIBRARIES)
        MESSAGE(STATUS "Found Python library: ${PYTHON_LIBRARIES}")
    ELSE()
        MESSAGE(STATUS "Could NOT find Python library")
    ENDIF()

    ######################################################################
    #
    #      find boost::python library
    #
    ######################################################################
    # 'FIND_PACKAGE(Boost COMPONENTS python)' is unreliable because it often selects
    # boost_python for the wrong Python version
    IF(Boost_FOUND)
        IF(Boost_USE_MULTITHREADED)
            # define names for thread-safe library variants
            SET(BOOST_PYTHON_NAMES
                    boost_python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}-mt
                    boost_python${PYTHON_VERSION_MAJOR}-mt
                    boost_python-mt)
        ENDIF()
        # define names for boost_python library variants
        # (may or may not be thread-safe)
        SET(BOOST_PYTHON_NAMES ${BOOST_PYTHON_NAMES}
                # Linux with multiple Python versions
                boost_python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}
                # Mac with Python 3
                boost_python${PYTHON_VERSION_MAJOR}
                # default
                boost_python)

        FIND_LIBRARY(Boost_PYTHON_LIBRARY
                     NAMES ${BOOST_PYTHON_NAMES}
                     HINTS "${Boost_LIBRARY_DIR}"
                     DOC "boost_python libraries")
    ENDIF()

    if(Boost_PYTHON_LIBRARY)
        MESSAGE(STATUS "Found boost_python library: ${Boost_PYTHON_LIBRARY}")
    else()
        MESSAGE(STATUS "Could NOT find boost_python library")
    endif()

    ######################################################################
    #
    #      find default install directory for Python modules
    #      (usually PYTHONDIR/Lib/site-packages)
    #
    ######################################################################
    IF(NOT DEFINED VIGRANUMPY_INSTALL_DIR OR VIGRANUMPY_INSTALL_DIR MATCHES "^$")
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                         "from distutils.sysconfig import *; print(get_python_lib(1))"
                          OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
        FILE(TO_CMAKE_PATH ${PYTHON_SITE_PACKAGES} VIGRANUMPY_INSTALL_DIR)
    ENDIF()
    SET(VIGRANUMPY_INSTALL_DIR ${VIGRANUMPY_INSTALL_DIR}
        CACHE PATH "where to install the VIGRA Python package" FORCE)
    # this is the install path relative to CMAKE_INSTALL_PREFIX,
    # use this in INSTALL() commands to get packaging right
    FILE(RELATIVE_PATH VIGRANUMPY_INSTALL_DIR ${CMAKE_INSTALL_PREFIX} ${VIGRANUMPY_INSTALL_DIR})

    ######################################################################
    #
    #      find numpy include directory
    #      (usually below PYTHONDIR/Lib/site-packages/numpy)
    #
    ######################################################################
    IF(NOT PYTHON_NUMPY_INCLUDE_DIR)
        # Note: we must suppress possible output of the 'from numpy... import *' command,
        #       because the output cannot be interpreted correctly otherwise
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                         "import sys, os; sys.stdout = open(os.devnull, 'w'); from numpy.distutils.misc_util import *; sys.__stdout__.write(' '.join(get_numpy_include_dirs()))"
                          RESULT_VARIABLE PYTHON_NUMPY_NOT_FOUND
                          OUTPUT_VARIABLE PYTHON_NUMPY_INCLUDE_DIR
                          OUTPUT_STRIP_TRAILING_WHITESPACE)
        IF(NOT PYTHON_NUMPY_NOT_FOUND)
            FILE(TO_CMAKE_PATH ${PYTHON_NUMPY_INCLUDE_DIR} PYTHON_NUMPY_INCLUDE_DIR)
        ELSE()
            SET(PYTHON_NUMPY_INCLUDE_DIR "PYTHON_NUMPY_INCLUDE_DIR-NOTFOUND")
        ENDIF()
    ENDIF()

    SET(PYTHON_NUMPY_INCLUDE_DIR ${PYTHON_NUMPY_INCLUDE_DIR}
        CACHE PATH "Path to numpy include files" FORCE)
    IF(PYTHON_NUMPY_INCLUDE_DIR)
        MESSAGE(STATUS "Searching for Python numpy: ok")
    ELSE()
        MESSAGE(STATUS "Could NOT find Python numpy ('import numpy.distutils.misc_util' failed)")
    ENDIF()

    ######################################################################
    #
    #      check if nosetests are installed
    #
    ######################################################################
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
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
    find_program(PYTHON_SPHINX sphinx-build "${PYTHON_PREFIX}/Scripts")

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
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                     "import sys; p = sys.platform; print('windows') if p.startswith('win') else p"
                      OUTPUT_VARIABLE PYTHON_PLATFORM OUTPUT_STRIP_TRAILING_WHITESPACE)

    ######################################################################
    #
    #      set outputs
    #
    ######################################################################
    INCLUDE(FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(VIGRANUMPY_DEPENDENCIES DEFAULT_MSG
                         PYTHONINTERP_FOUND PYTHON_INCLUDE_PATH PYTHON_LIBRARIES
                         Boost_PYTHON_LIBRARY PYTHON_NUMPY_INCLUDE_DIR VIGRANUMPY_INSTALL_DIR)

    IF(NOT VIGRANUMPY_INCLUDE_DIRS OR VIGRANUMPY_INCLUDE_DIRS MATCHES "-NOTFOUND")
        #note that the numpy include dir is set _before_ the python include dir, such that
        #installing a more recent version of numpy on top of an existing python installation
        #works (otherwise, numpy includes are picked up from ${PYTHON_INCLUDE_PATH}/numpy )
        SET(VIGRANUMPY_INCLUDE_DIRS ${PYTHON_NUMPY_INCLUDE_DIR} ${PYTHON_INCLUDE_PATH} ${Boost_INCLUDE_DIR})
    ENDIF()
    SET(VIGRANUMPY_INCLUDE_DIRS ${VIGRANUMPY_INCLUDE_DIRS}
        CACHE PATH "include directories needed by VIGRA Python bindings"
        FORCE)
    SET(VIGRANUMPY_LIBRARIES ${PYTHON_LIBRARIES} ${Boost_PYTHON_LIBRARY})
ELSE()
    MESSAGE(STATUS "Python not found. Make sure that Python is in your PATH or use 'cmake-gui' to set the PYTHON_EXECUTABLE variable manually.")
ENDIF()
