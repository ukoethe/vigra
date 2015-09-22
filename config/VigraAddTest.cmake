# Enhanced version of ADD_TEST
#
# Usage of this module:
#
#     VIGRA_ADD_TEST(target [[SOURCES] source1 source2 ...] [LIBRARIES lib1 lib2 ...])
#     VIGRA_COPY_TEST_DATA(datafile1 datafile2 ...)
#
# The function VIGRA_ADD_TEST
# * creates a new executable for 'target', using the given sources and libraries
# * makes the global target 'test' depend on the new 'target' (target 'test' must already exist)
# * installs a post-build event that runs the test automatically after linking
#
# The function VIGRA_COPY_TEST_DATA copies the given files from the current source directory
# to the corresponding binary directory.
#
MACRO(VIGRA_NATIVE_PATH out in)
    if(CMAKE_MAJOR_VERSION STREQUAL "3")
        cmake_policy(SET CMP0026 OLD)
    endif()
    file(TO_CMAKE_PATH "${in}" ${out})
    IF(NOT CMAKE_CFG_INTDIR STREQUAL ".")
        STRING(REGEX REPLACE "\\$\\([^\\)]*\\)" "%CONFIGURATION%" ${out} "${${out}}")
    ENDIF()
    IF(MINGW)
        # turn "c:/" into "/c/"
        STRING(REGEX REPLACE "^([a-zA-Z]):" "/\\1" ${out} "${${out}}")
    ENDIF()
    file(TO_NATIVE_PATH "${${out}}" ${out})
ENDMACRO(VIGRA_NATIVE_PATH)

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
    
    FILE(GLOB TESTSUCCESS_FOUND ${CMAKE_CURRENT_BINARY_DIR}/testsuccess.cxx)
    IF(NOT TESTSUCCESS_FOUND)
        FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/testsuccess.cxx 
         "// auto-generated dummy file to force re-execution of failed tests
")
    ENDIF()
    SET(SOURCES ${SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/testsuccess.cxx)
    
    # configure the target 
    IF(AUTOBUILD_TESTS)
        ADD_EXECUTABLE(${target} ${SOURCES})
    ELSE()
        ADD_EXECUTABLE(${target} EXCLUDE_FROM_ALL ${SOURCES})
    ENDIF()
    
    ADD_DEPENDENCIES(check_cpp ${target})
    ADD_DEPENDENCIES(ctest ${target})
    if(DEFINED LIBRARIES)
        TARGET_LINK_LIBRARIES(${target} ${LIBRARIES})
    endif()
    
    # find the test executable
    GET_TARGET_PROPERTY(${target}_executable ${target} LOCATION)
    VIGRA_NATIVE_PATH(VIGRA_TEST_EXECUTABLE ${${target}_executable})
    
    # Windows: set the DLL path
    set(VIGRA_PATH "")
    IF(MSVC)
        SET(PATHSEP ";")
    ELSE()
        SET(PATHSEP ":")
    ENDIF()
    FOREACH(lib ${LIBRARIES})
        GET_TARGET_PROPERTY(p ${lib} LOCATION)
        if(p)
            GET_FILENAME_COMPONENT(p ${p} PATH)
            VIGRA_NATIVE_PATH(p ${p})
            SET(VIGRA_PATH  "${p}${PATHSEP}${VIGRA_PATH}")
        endif()
    ENDFOREACH(lib)
    
    VIGRA_NATIVE_PATH(VIGRA_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
    IF(MSVC)
        SET(VIGRA_RUN_TEST "${CMAKE_CURRENT_BINARY_DIR}/run_${target}.bat")
        SET(VIGRA_TEST_EXECUTABLE "\"${VIGRA_TEST_EXECUTABLE}\"")  # take care of paths with spaces
        CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/config/run_test.bat.in
                       ${VIGRA_RUN_TEST}
                       @ONLY)
    ELSE()
        IF(VIGRA_RUN_TESTS_DIRECTLY)
          SET(VIGRA_RUN_TEST "${CMAKE_CURRENT_BINARY_DIR}/${target}")
        ELSE()
          SET(VIGRA_RUN_TEST "${CMAKE_CURRENT_BINARY_DIR}/run_${target}.sh")
          CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/config/run_test.sh.in
                         ${VIGRA_RUN_TEST}
                         @ONLY)
          EXECUTE_PROCESS(COMMAND chmod u+x ${VIGRA_RUN_TEST} OUTPUT_QUIET ERROR_QUIET)
        ENDIF()
    ENDIF()
    
    # register the test execution command
    IF(NOT CMAKE_CFG_INTDIR STREQUAL ".")
        SET(VIGRA_CONFIGURATION ${CMAKE_CFG_INTDIR})
    ELSE()
        SET(VIGRA_CONFIGURATION)
    ENDIF()
    
    IF(AUTOEXEC_TESTS)
        add_custom_command(
            TARGET ${target}
            POST_BUILD
            COMMAND ${VIGRA_RUN_TEST} ARGS ${VIGRA_CONFIGURATION}
            COMMENT "Running ${target}")
    ENDIF()
    
    ADD_TEST(${target} ${VIGRA_RUN_TEST} ${VIGRA_CONFIGURATION})
    
    IF(WITH_VALGRIND AND VALGRIND_EXECUTABLE)
        IF(VALGRIND_SUPPRESSION_FILE)
            SET(VALGRIND_SUPPRESSION "--suppressions=${VALGRIND_SUPPRESSION_FILE}")
        ELSE()
            SET(VALGRIND_SUPPRESSION)
        ENDIF()
        ADD_TEST(${target}_valgrind 
                ${VALGRIND_EXECUTABLE} 
                ${VALGRIND_SUPPRESSION} 
                --error-exitcode=1 
                ${${target}_executable})
    ENDIF()

ENDFUNCTION(VIGRA_ADD_TEST)

MACRO(VIGRA_COPY_TEST_DATA)
    FOREACH(test_data ${ARGN})
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${test_data} 
                       ${CMAKE_CURRENT_BINARY_DIR}/${test_data}
                       COPYONLY)
    ENDFOREACH(test_data)
ENDMACRO(VIGRA_COPY_TEST_DATA)
