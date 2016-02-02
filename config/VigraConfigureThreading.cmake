macro(VIGRA_CONFIGURE_THREADING)
    get_property(THREADING_IMPLEMENTATION GLOBAL PROPERTY THREADING_IMPLEMENTATION)
    if(NOT THREADING_IMPLEMENTATION)
        message(FATAL_ERROR "Threading implementation not detected yet: "
                            "Your CMakeLists file should not call VIGRA_CONFIGURE_THREADING() until after VigraDetectThreading "
                            "has been included from the root CMakeLists.txt file.")
    endif()

    set(THREADING_FOUND 1)
    if(THREADING_IMPLEMENTATION MATCHES "std")
        # Great, we can use the std library.  Nothing to do here...
    elseif(THREADING_IMPLEMENTATION MATCHES "none")
        if("${ARGN}" MATCHES "REQUIRED")
            message(FATAL_ERROR
                "***********************************************\n"
                "${THREADING_TEST_OUTPUT}\n"
                "***********************************************\n"
                "C++ threading is required, but your compiler doesn't support it. \n"
                "For reference, the failed compiler output is shown above.")
         else()
            add_definitions(-DVIGRA_SINGLE_THREADED)
            set(THREADING_FOUND 0)
        endif()
    else()
        message(FATAL_ERROR "CMakeLists bug: Unknown threading implementation")
    endif()
endmacro(VIGRA_CONFIGURE_THREADING)
