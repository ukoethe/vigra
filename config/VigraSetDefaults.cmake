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
