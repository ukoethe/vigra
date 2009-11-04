IF(NOT DEFINED WITH_VIGRANUMPY)
    SET(WITH_VIGRANUMPY "ON")
ENDIF()
SET(WITH_VIGRANUMPY ${WITH_VIGRANUMPY}
    CACHE BOOL "Build VIGRA Python bindings ?"
    FORCE)
    
IF(NOT DEFINED WITH_VALGRIND)
    SET(WITH_VALGRIND "OFF")
ENDIF()
SET(WITH_VALGRIND ${WITH_VALGRIND}
    CACHE BOOL "Perform valgrind memory testing upon 'make ctest' ?"
    FORCE)
    
IF(NOT DEFINED DEPENDENCY_SEARCH_PREFIX)
    SET(DEPENDENCY_SEARCH_PREFIX "")
ENDIF()    
SET(DEPENDENCY_SEARCH_PREFIX ${DEPENDENCY_SEARCH_PREFIX}
    CACHE PATH "Additional search prefixes (used by Find... macros)."
    FORCE)

IF(NOT DEFINED AUTOEXEC_TESTS)
    SET(AUTOEXEC_TESTS "ON")
ENDIF()    
SET(AUTOEXEC_TESTS ${AUTOEXEC_TESTS}
    CACHE BOOL "Automatically execute each test after compilation ?"
    FORCE)

IF(NOT DEFINED AUTOBUILD_TESTS)
    SET(AUTOBUILD_TESTS "OFF")
ENDIF()    
SET(AUTOBUILD_TESTS ${AUTOBUILD_TESTS}
    CACHE BOOL "Compile tests as part of target 'all' (resp. 'ALL_BUILD') ?"
    FORCE)

IF(NOT VIGRA_DEFAULTS_INIT)
    SET(VIGRA_DEFAULTS_INIT TRUE CACHE INTERNAL "initial flags set")

    IF (NOT MSVC AND NOT CMAKE_BUILD_TYPE)
        SET(CMAKE_BUILD_TYPE "Release"
            CACHE STRING "Choose the type of build, options are None Release Debug RelWithDebInfo MinSizeRel." FORCE)
    ENDIF ()
    
#    # initial compiler flags can be set here, this is only
#    # executed once in the first configure run.
#    IF(CMAKE_COMPILER_IS_GNUCXX)
#        IF(NOT CMAKE_C_FLAGS)
#            SET(CMAKE_C_FLAGS
#                "-W -Wall"
#                CACHE STRING
#                "Flags used by the compiler during all build types."
#                FORCE
#            )
#        ENDIF(NOT CMAKE_C_FLAGS)
#        IF(NOT CMAKE_CXX_FLAGS)
#            SET(CMAKE_CXX_FLAGS
#                "-W -Wall"
#                CACHE STRING
#                "Flags used by the compiler during all build types."
#                FORCE
#            )
#        ENDIF(NOT CMAKE_CXX_FLAGS)
#    ENDIF(CMAKE_COMPILER_IS_GNUCXX)
ENDIF(NOT VIGRA_DEFAULTS_INIT)

