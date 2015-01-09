macro(VIGRA_CONFIGURE_THREADING)
    get_property(THREADING_IMPLEMENTATION GLOBAL PROPERTY THREADING_IMPLEMENTATION)
    if(NOT THREADING_IMPLEMENTATION)
        message(FATAL_ERROR "Threading implementation not detected yet: "
                            "Your CMakeLists file should not call VIGRA_CONFIGURE_THREADING() until after VigraDetectThreading " 
                            "has been included from the root CMakeLists.txt file.")
    endif()

    set(THREADING_FOUND 1)
    if(THREADING_IMPLEMENTATION MATCHES "boost")
        include_directories(${Boost_INCLUDE_DIR})
        set(THREADING_LIBRARIES ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${Boost_CHRONO_LIBRARY})
    elseif(THREADING_IMPLEMENTATION MATCHES "std")
        # Great, we can use the std library.  Nothing to do here...
    elseif(THREADING_IMPLEMENTATION MATCHES "none")
        if("${ARGN}" MATCHES "REQUIRED")
            message(FATAL_ERROR 
                "***********************************************\n"
                "${THREADING_TEST_OUTPUT}\n"
                "***********************************************\n"
                "C++ threading is required, but your compiler doesn't support it. \n"
                "For reference, the failed compiler output is shown above. \n"
                "Try to run cmake with '-DWITH_BOOST_THREAD=1'.\n")
        else()
            message("(Compiling single threaded, consider running 'cmake -DWITH_BOOST_THREAD=1')")
            add_definitions(-DVIGRA_SINGLE_THREADED)
            set(THREADING_FOUND 0)
        endif()
    else()
        message(FATAL_ERROR "CMakeLists bug: Unknown threading implementation")
    endif()
endmacro(VIGRA_CONFIGURE_THREADING)
