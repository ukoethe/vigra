# - Adds the right OpenMP flags to the CXX compiler
# The following variables are set:
#   OpenMP_CXX_FLAGS - flags to add to the CXX compiler for OpenMP support
#   OPENMP_FOUND - true if openmp is detected

FIND_PACKAGE(OpenMP)

# handle OPENMP flags
IF(OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()