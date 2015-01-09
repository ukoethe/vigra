# Check for a valid isnan(int) implementation
# (Some old versions of libc++ don't provide isnan() for integer types.)
try_compile(ISNAN_INTEGERS_OKAY
    ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/config/check_std_isnan.cxx )

if (NOT ISNAN_INTEGERS_OKAY)
    ADD_DEFINITIONS(-DVIGRA_NO_INT_ISNAN)
endif()