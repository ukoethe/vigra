# - Find VIGRANUMPY_DEPENDENCIES
#
MESSAGE(STATUS "Checking VIGRANUMPY_DEPENDENCIES")

IF(NOT Python_Interpreter_FOUND)
    FIND_PACKAGE(Python COMPONENTS Interpreter Development NumPy)
ENDIF()

IF(Python_Interpreter_FOUND)
    MESSAGE(STATUS "Using Python ${Python_VERSION} at ${Python_EXECUTABLE}")
    MESSAGE(STATUS "Python_LIBRARIES ${Python_LIBRARIES}")

    ######################################################################
    #
    #      find Python prefix
    #
    ######################################################################
    execute_process(COMMAND ${Python_EXECUTABLE} -c
                     "import sys; print(sys.exec_prefix)"
                      OUTPUT_VARIABLE PYTHON_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)

    IF(Python_INCLUDE_DIRS)
        MESSAGE(STATUS "Found Python includes:  ${Python_INCLUDE_DIRS}")
    ELSE()
        MESSAGE(STATUS "Could NOT find Python includes")
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
                    boost_python-py${Python_VERSION_MAJOR}${Python_VERSION_MINOR}-mt
                    boost_python-${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}-mt
                    boost_python${Python_VERSION_MAJOR}-mt
                    boost_python-mt)
        ENDIF()
        IF(Boost_LIB_SUFFIX)
            SET(BOOST_PYTHON_NAMES ${BOOST_PYTHON_NAMES}
                # Windows with mangled library names
                boost_python${Python_VERSION_MAJOR}${Boost_LIB_SUFFIX}
                boost_python${Boost_LIB_SUFFIX})
        ENDIF()

        # define names for boost_python library variants
        # (may or may not be thread-safe)
        SET(BOOST_PYTHON_NAMES ${BOOST_PYTHON_NAMES}
                # Linux with multiple Python versions
                boost_python-py${Python_VERSION_MAJOR}${Python_VERSION_MINOR}
                # Gentoo
                boost_python-${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}
                # Mac with Python 3
                boost_python${Python_VERSION_MAJOR}
                # conda-forge
                boost_python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}
                # default
                boost_python)

        FIND_LIBRARY(Boost_PYTHON_LIBRARY
                     NAMES ${BOOST_PYTHON_NAMES}
                     NAMES_PER_DIR
                     HINTS "${Boost_LIBRARY_DIR}"
                     DOC "boost_python libraries")
    ENDIF()

    if(Boost_PYTHON_LIBRARY)
        MESSAGE(STATUS "Found boost_python library: ${Boost_PYTHON_LIBRARY}")
    else()
        MESSAGE(STATUS "Could NOT find boost_python library (looking for version ${Python_VERSION_MAJOR}.${Python_VERSION_MINOR})")
    endif()

    ######################################################################
    #
    #      find default install directory for Python modules
    #      (usually PYTHONDIR/Lib/site-packages)
    #
    ######################################################################
    IF(NOT DEFINED VIGRANUMPY_INSTALL_DIR OR VIGRANUMPY_INSTALL_DIR MATCHES "^$")
        execute_process(COMMAND ${Python_EXECUTABLE} -c
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
    #      check if nosetests are installed
    #
    ######################################################################
    execute_process(COMMAND ${Python_EXECUTABLE} -c
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
    execute_process(COMMAND ${Python_EXECUTABLE} -c
                     "import sys; p = sys.platform; print('windows') if p.startswith('win') else p"
                      OUTPUT_VARIABLE PYTHON_PLATFORM OUTPUT_STRIP_TRAILING_WHITESPACE)

    ######################################################################
    #
    #      set outputs
    #
    ######################################################################
    INCLUDE(FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(VIGRANUMPY_DEPENDENCIES DEFAULT_MSG
                         Python_Interpreter_FOUND Python_INCLUDE_DIRS Python_LIBRARIES
                         Boost_PYTHON_LIBRARY Python_NumPy_INCLUDE_DIRS VIGRANUMPY_INSTALL_DIR)

    IF(NOT VIGRANUMPY_INCLUDE_DIRS OR VIGRANUMPY_INCLUDE_DIRS MATCHES "-NOTFOUND")
        #note that the numpy include dir is set _before_ the python include dir, such that
        #installing a more recent version of numpy on top of an existing python installation
        #works (otherwise, numpy includes are picked up from ${PYTHON_INCLUDE_PATH}/numpy )
        SET(VIGRANUMPY_INCLUDE_DIRS ${Python_NumPy_INCLUDE_DIRS} ${Python_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
    ENDIF()
    SET(VIGRANUMPY_INCLUDE_DIRS ${VIGRANUMPY_INCLUDE_DIRS}
        CACHE PATH "include directories needed by VIGRA Python bindings"
        FORCE)
    SET(VIGRANUMPY_LIBRARIES ${Python_LIBRARIES} ${Boost_PYTHON_LIBRARY})

    if(TEST_VIGRANUMPY AND NOT VIGRANUMPY_DEPENDENCIES_FOUND)
        MESSAGE(FATAL_ERROR "  vigranumpy dependencies NOT found while TEST_VIGRANUMPY=1")
    endif()
    if(TEST_VIGRANUMPY AND PYTHON_NOSETESTS_NOT_FOUND)
        MESSAGE(FATAL_ERROR "  nosetests NOT found while TEST_VIGRANUMPY=1")
    endif()
ELSE()
    MESSAGE(STATUS "Python not found. Make sure that Python is in your PATH or set the appropriate variables as described in https://cmake.org/cmake/help/latest/module/FindPython.html")
ENDIF()
