MACRO(VIGRA_CONFIGURE_THREADING)
    # We store a global property to avoid running this macro repeatedly. 
    # (Avoid running the compiler repeatedly.)
    get_property(THREADING_ALREADY_TESTED GLOBAL PROPERTY THREADING_FOUND SET)
    if (NOT THREADING_ALREADY_TESTED)
        if(WITH_BOOST_THREAD)
            if(Boost_FOUND)
                # See boost library list, below.
                MESSAGE(STATUS "Checking for threading support:   boost::thread")
                set_property(GLOBAL PROPERTY THREADING_FOUND 1)
            else()
                MESSAGE(FATAL_ERROR "Cmake was run with '-DWITH_BOOST_THREAD=1', but required boost components were not found.")
            endif()
        else()
            TRY_COMPILE(STD_THREADING_FOUND 
                ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/config/check_std_thread.cxx
                OUTPUT_VARIABLE THREADING_TEST_OUTPUT )
            
            set_property(GLOBAL PROPERTY THREADING_TEST_OUTPUT "${THREADING_TEST_OUTPUT}")
    
            if (STD_THREADING_FOUND)
                MESSAGE(STATUS "Checking for threading support:   std::thread")
                set_property(GLOBAL PROPERTY THREADING_FOUND 1)

                if (CMAKE_COMPILER_IS_GNUCXX AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.7.0")
                    SET(CMAKE_CXX_FLAGS "-pthread -std=c++0x ${CMAKE_CXX_FLAGS}")
                elseif(CMAKE_COMPILER_IS_GNUCXX)
                    SET(CMAKE_CXX_FLAGS "-pthread -std=c++11 ${CMAKE_CXX_FLAGS}")
                elseif(NOT MSVC)
                    SET(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
                endif()
            else()
                MESSAGE(STATUS "Checking for threading support:   NOT FOUND")
                set_property(GLOBAL PROPERTY THREADING_FOUND 0)

                if(NOT "${ARGN}" MATCHES "REQUIRED")
                    MESSAGE("(Compiling single threaded, consider running 'cmake -DWITH_BOOST_THREAD=1')")
                    ADD_DEFINITIONS(-DVIGRA_SINGLE_THREADED)
                endif()
            endif()
        endif()
    endif()

    # If at any module calls this macro with "REQUIRED" as the argument, 
    #  we fail unless a threading implementation has been found.
    get_property(THREADING_FOUND GLOBAL PROPERTY THREADING_FOUND)
    if((NOT THREADING_FOUND) AND ("${ARGN}" MATCHES "REQUIRED"))
        get_property(THREADING_TEST_OUTPUT GLOBAL PROPERTY THREADING_TEST_OUTPUT)
        MESSAGE(FATAL_ERROR 
            "***********************************************\n"
            "${THREADING_TEST_OUTPUT}\n"
            "***********************************************\n"
            "C++ threading is required, but your compiler doesn't support it. \n"
            "For reference, the failed compiler output is shown above. \n"
            "Try to run cmake with '-DWITH_BOOST_THREAD=1'.\n")
    endif()    

    # We list the boost libraries down here, 
    #  so they are available in the caller's local scope
    if(WITH_BOOST_THREAD)
        INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
        SET(THREADING_LIBRARIES ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${Boost_CHRONO_LIBRARY})
    endif()

ENDMACRO(VIGRA_CONFIGURE_THREADING)
