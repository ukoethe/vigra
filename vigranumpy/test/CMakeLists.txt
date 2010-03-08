INCLUDE_DIRECTORIES(${VIGRANUMPY_INCLUDE_DIRS})

FILE(GLOB TESTSUCCESS_FOUND ${CMAKE_CURRENT_BINARY_DIR}/testsuccess.cxx)
IF(NOT TESTSUCCESS_FOUND)
    FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/testsuccess.cxx 
     "// auto-generated dummy file to force re-execution of failed tests
")
ENDIF()

ADD_LIBRARY(vigranumpytest SHARED EXCLUDE_FROM_ALL
             vigranumpytest.cxx ${CMAKE_CURRENT_BINARY_DIR}/testsuccess.cxx)

IF(PYTHON_PLATFORM MATCHES "^win32$")
    SET_TARGET_PROPERTIES(vigranumpytest PROPERTIES PREFIX "" SUFFIX  ".pyd")
ELSE()
    SET_TARGET_PROPERTIES(vigranumpytest PROPERTIES PREFIX "")
ENDIF()

TARGET_LINK_LIBRARIES(vigranumpytest ${VIGRANUMPY_LIBRARIES} vigranumpy_core)
ADD_DEPENDENCIES(check_python vigranumpytest)
ADD_DEPENDENCIES(ctest vigranumpytest)
ADD_DEPENDENCIES(vigranumpytest vigranumpy)

SET(TEST_SCRIPTS
    test_arraytypes.py
    test_impex.py
    test1.py
    test2.py
    test3.py)

INCLUDE(VigraAddTest)
VIGRA_COPY_TEST_DATA(${TEST_SCRIPTS})

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

 
        VIGRA_NATIVE_PATH(PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE})
        
        # Windows: set the DLL path
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
