# - Adds the right OpenMP flags to the CXX compiler

FIND_PACKAGE(OpenMP)

IF(OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()