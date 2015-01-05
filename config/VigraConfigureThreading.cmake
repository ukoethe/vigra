MACRO(VIGRA_CONFIGURE_THREADING)
    set(THREADING_FOUND 0)
    if(((APPLE AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.7.0") OR MINGW)
         AND NOT WITH_BOOST_THREAD)
        TRY_COMPILE(THREADING_FOUND 
            ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/config/check_std_thread.cxx
            OUTPUT_VARIABLE THREADING_TEST_OUTPUT )
        if (NOT THREADING_FOUND)
            if(${ARGN} MATCHES "REQUIRED")
                MESSAGE(FATAL_ERROR "Compiler does not support C++ threading. Try to run cmake with '-DWITH_BOOST_THREAD=1'.\n"
                "For reference, the failed compiler output is shown below."
                "${THREADING_TEST_OUTPUT}\n\n"
                "Compiler does not support C++ threading. Try to run cmake with '-DWITH_BOOST_THREAD=1'.\n"
                "For reference, the failed compiler output is shown above.\n"
                )
            else()
                ADD_DEFINITIONS(-DVIGRA_SINGLE_THREADED)
            endif()
        endif()
    else()
        if(WITH_BOOST_THREAD)
            if(Boost_FOUND)
                set(THREADING_FOUND 1)
                INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
                SET(THREADING_LIBRARIES ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${Boost_CHRONO_LIBRARY})
            elseif(${ARGV1} MATCHES "REQUIRED")  
                MESSAGE(FATAL_ERROR "Cmake was run with '-DWITH_BOOST_THREAD=1', but required boost components were not found.")
            else()
                ADD_DEFINITIONS(-DVIGRA_SINGLE_THREADED)
            endif()
        else()
            # note: we only support MSVC versions that have threading
            set(THREADING_FOUND 1)
            IF (CMAKE_COMPILER_IS_GNUCXX AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.7.0")
                SET(CMAKE_CXX_FLAGS "-pthread -std=c++0x ${CMAKE_CXX_FLAGS}")
            elseif(CMAKE_COMPILER_IS_GNUCXX)
                SET(CMAKE_CXX_FLAGS "-pthread -std=c++11 ${CMAKE_CXX_FLAGS}")
            elseif(NOT MSVC)
                SET(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
            endif()
        endif()
    endif()
ENDMACRO(VIGRA_CONFIGURE_THREADING)
