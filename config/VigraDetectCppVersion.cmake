macro(VIGRA_DETECT_CPP_VERSION)

try_run(RUN_RESULT COMPILE_SUCCEEDED
        ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/config/output_cplusplus_version.cxx
        COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT_FROM_CPP_DETECT
        RUN_OUTPUT_VARIABLE VIGRA_CPP_VERSION)

    if(NOT COMPILE_SUCCEEDED)
        message(FATAL_ERROR "Failed to detect c++ version with a simple test program!  Compiler Output is below.\n"
                "${COMPILE_OUTPUT_FROM_CPP_DETECT}")
    endif()

    if(RUN_RESULT)
        message(FATAL_ERROR "Failed to detect c++ version with a simple test program!\n"
			    "Test program compiled, but did not execute cleanly. Run output is shown below.\n"
 			    "${VIGRA_CPP_VERSION}")
    endif()

    message(STATUS "Detected C++ version: ${VIGRA_CPP_VERSION}")

endmacro(VIGRA_DETECT_CPP_VERSION)
