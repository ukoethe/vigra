INCLUDE_DIRECTORIES(${VIGRANUMPY_INCLUDE_DIR})

ADD_LIBRARY(vigranumpytest SHARED EXCLUDE_FROM_ALL
            vigranumpytest.cxx)

IF(PYTHON_PLATFORM MATCHES "^win32$")
    SET_TARGET_PROPERTIES(vigranumpytest PROPERTIES PREFIX "" SUFFIX  ".pyd")
ELSE()
    SET_TARGET_PROPERTIES(vigranumpytest PROPERTIES PREFIX "")
ENDIF()

TARGET_LINK_LIBRARIES(vigranumpytest ${VIGRANUMPY_LIBRARIES} vigranumpycmodule)
ADD_DEPENDENCIES(check vigranumpytest)

SET(TEST_SCRIPTS
    test_arraytypes.py
    test_impex.py
    test_christoph.py
    test_nathan.py
    test_stephan.py)

INCLUDE(VigraAddTest)
VIGRA_COPY_TEST_DATA(${TEST_SCRIPTS})

MACRO(DEPENDENCY_PATH variable target)
    GET_TARGET_PROPERTY(${variable} ${target} LOCATION)
    STRING(REGEX REPLACE "(/\\$\\(OutDir\\)/[^/]*|/[^/]*)$" "" ${variable} ${${variable}}) # get path prefix
    FILE(TO_NATIVE_PATH ${${variable}} ${variable})
ENDMACRO(DEPENDENCY_PATH)

DEPENDENCY_PATH(VIGRAIMPEX_PATH vigraimpex)
DEPENDENCY_PATH(VIGRANUMPYCMODULE_PATH vigranumpycmodule)
DEPENDENCY_PATH(VIGRANUMPYTEST_PATH vigranumpytest)

FILE(TO_CMAKE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../lib VIGRANUMPYSCRIPTS_PATH)
FILE(TO_NATIVE_PATH ${VIGRANUMPYSCRIPTS_PATH} VIGRANUMPYSCRIPTS_PATH)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/set_paths.py 
               ${CMAKE_CURRENT_BINARY_DIR}/set_paths.py 
               @ONLY)

IF(AUTOEXEC_TESTS)
    IF(NOSETESTS_EXECUTABLE)
        IF (MSVC)
            SET(TEST_ARGS . $(OutDir))
        ELSE()
            SET(TEST_ARGS .)
        ENDIF()

        add_custom_command(
            TARGET vigranumpytest
            POST_BUILD
            COMMAND ${NOSETESTS_EXECUTABLE} ARGS ${TEST_ARGS}
            COMMENT "Running vigranumpy tests")
    ELSE()
        MESSAGE(STATUS "  vigranumpy tests will NOT be executed (nosetests missing)")
    ENDIF()
ENDIF()
