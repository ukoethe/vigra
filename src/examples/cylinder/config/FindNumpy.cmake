#adapted from vigra source code

IF(NOT PYTHON_NUMPY_INCLUDE_DIR)
    execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c 
                        "from numpy.distutils.misc_util import *; print ' '.join(get_numpy_include_dirs())" 
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
    SET(NUMPY_FOUND TRUE)

    IF (NOT Numpy_FIND_QUIETLY)
      MESSAGE(STATUS "Found Numpy")
      MESSAGE(STATUS "  > includes:      ${PYTHON_NUMPY_INCLUDE_DIR}")
    ENDIF()
ELSE()
    IF(Numpy_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find Numpy")
    ENDIF()
ENDIF()
