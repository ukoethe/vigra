# Enhanced version of ADD_TEST
#
# Usage of this module:
#
#     VIGRA_ADD_TEST(target [[SOURCES] source1 source2 ...] [LIBRARIES lib1 lib2 ...])
#     VIGRA_COPY_TEST_DATA(datafile1 datafile2 ...)
#
# The function VIGRA_ADD_TEST
# * creates a new executable for 'target', using the given sources and libraries
# * makes the global target 'test' depend on the new 'target'
# * installs a post-build event that runs the test automatically after linking
#
# The function VIGRA_COPY_TEST_DATA copies the given files from the current source directory
# to the corresponding binary directory.
#
INCLUDE(testOrDelete RESULT_VARIABLE TEST_OR_DELETE)

IF (CMAKE_GENERATOR MATCHES "Visual Studio")
   STRING(REGEX REPLACE "\\.cmake$" ".bat" TEST_OR_DELETE ${TEST_OR_DELETE})
ELSE ()
   STRING(REGEX REPLACE "\\.cmake$" ".sh" TEST_OR_DELETE ${TEST_OR_DELETE})
ENDIF ()

if (WIN32 AND ${CMAKE_C_COMPILER} MATCHES gcc.exe)
    SET(IS_CYGWIN 1)
else()
    SET(IS_CYGWIN 0)
endif()

FUNCTION(VIGRA_ADD_TEST target)
    # parse the args
    set(v SOURCES)
    foreach(i ${ARGN})
        if(${i} MATCHES "^SOURCES$")
            set(v SOURCES)
        elseif(${i} MATCHES "^LIBRARIES$")
            set(v LIBRARIES)
        else()
            set(${v} ${${v}} ${i})
        endif()
    endforeach(i)
    
    # configure the target
    ADD_EXECUTABLE(${target} ${SOURCES})
    ADD_DEPENDENCIES(test ${target})
    if(DEFINED LIBRARIES)
        TARGET_LINK_LIBRARIES(${target} ${LIBRARIES})
    endif()
    
    # find the test executable
    GET_TARGET_PROPERTY(${target}_executable ${target} LOCATION)
    file(TO_NATIVE_PATH ${${target}_executable} ${target}_executable)
    
    # cygwin: set the DLL path (not needed for other platforms
    set(path "")
    if(IS_CYGWIN)
        FOREACH(lib ${LIBRARIES})
            GET_TARGET_PROPERTY(p ${lib} LOCATION)
            STRING(REGEX REPLACE "/[^/]*$" "" p ${p})
            if(NOT ${p} MATCHES "NOTFOUND")
                set(path "${path}${p}:")
            endif()
        ENDFOREACH(lib)
    endif()
    
    # register the test execution command
    add_custom_command(
        TARGET ${target}
        POST_BUILD
        COMMAND ${TEST_OR_DELETE} ARGS ${${target}_executable} ${path} 
        COMMENT "Running tests")
ENDFUNCTION(VIGRA_ADD_TEST)

MACRO(VIGRA_COPY_TEST_DATA)
    FOREACH(test_data ${ARGN})
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${test_data} 
                       ${CMAKE_CURRENT_BINARY_DIR}/${test_data}
                       COPYONLY)
    ENDFOREACH(test_data)
ENDMACRO(VIGRA_COPY_TEST_DATA)
