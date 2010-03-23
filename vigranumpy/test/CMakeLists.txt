INCLUDE_DIRECTORIES(${VIGRANUMPY_INCLUDE_DIRS})

SET(TEST_SCRIPTS
    test_arraytypes.py
    test_impex.py
    test1.py
    test2.py
    test3.py)

# setup the file 'testsuccess.cxx' which will become out-of-date when the
# last test failed or when the tests code changed, so that tests are re-executed
SET(TEST_SUCCESS_FILE "${CMAKE_CURRENT_BINARY_DIR}/testsuccess.cxx")

# find the vigranumpy modules that we want to test, and make
# 'testsuccess.cxx' depend on them
# TODO: Does cmake provide a standard way to find those dependencies? 
GET_TARGET_PROPERTY(VIGRANUMPY_DEPENDS vigranumpy VIGRA_DEPENDS)
SEPARATE_ARGUMENTS(VIGRANUMPY_DEPENDS)

add_custom_command(
    OUTPUT ${TEST_SUCCESS_FILE}
    DEPENDS ${TEST_SCRIPTS} ${VIGRANUMPY_DEPENDS}
    COMMAND ${CMAKE_COMMAND}
    ARGS -E touch ${TEST_SUCCESS_FILE})

FILE(GLOB TESTSUCCESS_FOUND ${TEST_SUCCESS_FILE})
IF(NOT TESTSUCCESS_FOUND)
    FILE(WRITE ${TEST_SUCCESS_FILE}
     "// auto-generated dummy file to force re-execution of failed tests
")
ENDIF()

# setup the test target
ADD_LIBRARY(vigranumpytest SHARED EXCLUDE_FROM_ALL
             vigranumpytest.cxx ${TEST_SUCCESS_FILE} ${TEST_SCRIPTS})

IF(PYTHON_PLATFORM MATCHES "^win32$")
    SET_TARGET_PROPERTIES(vigranumpytest PROPERTIES PREFIX "" SUFFIX  ".pyd")
ELSE()
    SET_TARGET_PROPERTIES(vigranumpytest PROPERTIES PREFIX "")
ENDIF()

TARGET_LINK_LIBRARIES(vigranumpytest ${VIGRANUMPY_LIBRARIES} vigranumpy_core)
ADD_DEPENDENCIES(check_python vigranumpytest)
ADD_DEPENDENCIES(ctest vigranumpytest)
ADD_DEPENDENCIES(vigranumpytest vigranumpy)

# setup the file 'set_path.py' to ensure that tests are run against the binaries
# which are currently being built, not against (possibly outdated) installed modules
MACRO(DEPENDENCY_PATH variable target)
    GET_TARGET_PROPERTY(${variable} ${target} LOCATION)
    STRING(REGEX REPLACE "(/\\$\\([^\\)]*\\)/[^/]*|/[^/]*)$" "" ${variable} ${${variable}}) # get path prefix
    FILE(TO_NATIVE_PATH ${${variable}} ${variable})
ENDMACRO(DEPENDENCY_PATH)

DEPENDENCY_PATH(VIGRAIMPEX_PATH vigraimpex)
DEPENDENCY_PATH(VIGRANUMPYTEST_PATH vigranumpytest)

STRING(REGEX REPLACE "/vigra$" "" VIGRANUMPY_TMP_PATH ${vigranumpy_tmp_dir})
FILE(TO_NATIVE_PATH "${VIGRANUMPY_TMP_PATH}" VIGRANUMPY_TMP_PATH)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/set_paths.py 
                ${CMAKE_CURRENT_BINARY_DIR}/set_paths.py 
                @ONLY)

IF(AUTOEXEC_TESTS)
    IF(NOT PYTHON_NOSETESTS_NOT_FOUND)
     
        # copy the individual test scripts
        FOREACH(test_script ${TEST_SCRIPTS})
            ADD_CUSTOM_COMMAND(
                TARGET vigranumpytest
                POST_BUILD
                COMMAND ${CMAKE_COMMAND}
                ARGS -E copy_if_different
                ${CMAKE_CURRENT_SOURCE_DIR}/${test_script} ${CMAKE_CURRENT_BINARY_DIR}/${test_script}
                )
        ENDFOREACH(test_script)

        VIGRA_NATIVE_PATH(PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE})
        
        # Windows: Set the DLL path to the binaries currently being built
        #          (so that binaries which are already installed are ignored)
        #          This is not necessary on Linux, because the correct path
        #          is hard-coded into the module by means of 'rpath'
        set(VIGRA_PATH "")
        if(WIN32)
            IF(CYGWIN)
                SET(PATHSEP ":")
            ELSE()
                SET(PATHSEP ";")
            ENDIF()
            FOREACH(lib ${LIBRARIES})
                GET_TARGET_PROPERTY(p ${lib} LOCATION)
                if(p)
                    GET_FILENAME_COMPONENT(p ${p} PATH)
                    VIGRA_NATIVE_PATH(p ${p})
                    SET(VIGRA_PATH  "${p}${PATHSEP}${VIGRA_PATH}")
               endif()
            ENDFOREACH(lib)
        endif()
        
        # Setup the test execution script
        VIGRA_NATIVE_PATH(VIGRA_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
        SET(VIGRA_TEST_EXECUTABLE "${PYTHON_EXECUTABLE} -c \"import nose; nose.main()\" .")        
        IF(MSVC OR MINGW)
            SET(VIGRA_TEST_EXECUTABLE "${VIGRA_TEST_EXECUTABLE} %CONFIGURATION%")
            SET(VIGRA_RUN_TEST "${CMAKE_CURRENT_BINARY_DIR}/run_vigramumpytest.bat")
            CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/config/run_test.bat.in
                            ${VIGRA_RUN_TEST}
                            @ONLY)
        ELSE()            
            SET(VIGRA_RUN_TEST "${CMAKE_CURRENT_BINARY_DIR}/run_vigramumpytest.sh")
            CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/config/run_test.sh.in
                            ${VIGRA_RUN_TEST}
                           @ONLY)
            EXECUTE_PROCESS(COMMAND chmod u+x ${VIGRA_RUN_TEST} OUTPUT_QUIET ERROR_QUIET)
        ENDIF()

        # register the test execution command
        # VIGRA_CONFIGURATION is 'release' or 'debug' on Windows, nothing on Linux
        IF(NOT CMAKE_CFG_INTDIR STREQUAL ".")
            SET(VIGRA_CONFIGURATION ${CMAKE_CFG_INTDIR})
        ELSE()
            SET(VIGRA_CONFIGURATION)
        ENDIF()

        add_custom_command(
            TARGET vigranumpytest
            POST_BUILD
            COMMAND ${VIGRA_RUN_TEST} ARGS ${VIGRA_CONFIGURATION}
            COMMENT "Running vigranumpy tests")

    ELSE()
        MESSAGE(STATUS "  vigranumpy tests will NOT be executed (nosetests missing)")
    ENDIF()
ENDIF()
