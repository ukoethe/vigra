#
# This file determines which threading implementation, if any should be used for compiling vigra source files.
# Source modules should use VIGRA_CONFIGURE_THREADING() from VigraConfigureThreading.cmake to set includes/preprocessor definitions.
#

cmake_minimum_required(VERSION 2.8)

get_property(THREADING_IMPLEMENTATION GLOBAL PROPERTY THREADING_IMPLEMENTATION)
if(THREADING_IMPLEMENTATION)
    message(FATAL_ERROR "VigraDetectThreading.cmake should be included exactly once, from the root CMakeLists.txt file.  There should be no need to include it more than once.")
endif()

if(WITH_BOOST_THREAD)
    if(Boost_FOUND)
        # See boost library list, below.
        MESSAGE(STATUS "Checking for threading support:   boost::thread")
        set_property(GLOBAL PROPERTY THREADING_IMPLEMENTATION "boost")
    else()
        MESSAGE(FATAL_ERROR "Cmake was run with '-DWITH_BOOST_THREAD=1', but required boost components were not found.")
    endif()
else()
    # Save original flags before we add new ones, in case there's no need for the new ones.
    set(ORIG_CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})

    if (CMAKE_COMPILER_IS_GNUCXX AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.7.0")
        SET(CMAKE_CXX_FLAGS "-pthread -std=c++0x ${CMAKE_CXX_FLAGS}")
    elseif(CMAKE_COMPILER_IS_GNUCXX)
        SET(CMAKE_CXX_FLAGS "-pthread -std=c++11 ${CMAKE_CXX_FLAGS}")
    elseif(NOT MSVC)
        SET(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
    endif()

    TRY_COMPILE(STD_THREADING_FOUND 
        ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/config/check_std_thread.cxx
        OUTPUT_VARIABLE THREADING_TEST_OUTPUT )
    
    if(STD_THREADING_FOUND)
        MESSAGE(STATUS "Checking for threading support:   std::thread")
        set_property(GLOBAL PROPERTY THREADING_IMPLEMENTATION "std")
    else()
        MESSAGE(STATUS "Checking for threading support:   NOT FOUND")
        set_property(GLOBAL PROPERTY THREADING_IMPLEMENTATION "none")
        # Revert to old CXX_FLAGS
        set(CMAKE_CXX_FLAGS ${ORIG_CMAKE_CXX_FLAGS})
    endif()
endif()
