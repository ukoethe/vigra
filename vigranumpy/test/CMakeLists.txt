INCLUDE_DIRECTORIES(${VIGRANUMPY_INCLUDE_DIRS})

ADD_LIBRARY(vigranumpytest SHARED EXCLUDE_FROM_ALL
            vigranumpytest.cxx)

IF(PYTHON_PLATFORM MATCHES "^win32$")
    SET_TARGET_PROPERTIES(vigranumpytest PROPERTIES PREFIX "" SUFFIX  ".pyd")
ELSE()
    SET_TARGET_PROPERTIES(vigranumpytest PROPERTIES PREFIX "")
ENDIF()

TARGET_LINK_LIBRARIES(vigranumpytest ${VIGRANUMPY_LIBRARIES} vigranumpycore)
ADD_DEPENDENCIES(check vigranumpytest)

SET(TEST_SCRIPTS
    test_arraytypes.py
    test_impex.py
    test1.py
    test2.py
    test3.py)

INCLUDE(VigraAddTest)
VIGRA_COPY_TEST_DATA(${TEST_SCRIPTS})

MACRO(VIGRA_COPY_TEST_DATA)
    FOREACH(test_data ${ARGN})
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${test_data} 
                       ${CMAKE_CURRENT_BINARY_DIR}/${test_data}
                       COPYONLY)
    ENDFOREACH(test_data)
ENDMACRO(VIGRA_COPY_TEST_DATA)

#copy rest file into build-vigra module
SET(LIB_FILES __init__.py arraytypes.py ufunc.py)
FOREACH(lib_file ${LIB_FILES})
    ADD_CUSTOM_COMMAND(
        TARGET vigranumpytest
        PRE_BUILD
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy_if_different ${PROJECT_SOURCE_DIR}/vigranumpy/lib/${lib_file} ../src/vigra/${lib_file}
        )
ENDFOREACH(lib_file)

#copy vigranumpy files into build-vigra module
SET(MODULE_TARGETS vigranumpycore convolution segmentation edgedetection
    morphology classification tensor impex filters noise)

ADD_DEPENDENCIES(vigranumpytest vigranumpy)
FOREACH(module ${MODULE_TARGETS})
    GET_TARGET_PROPERTY(loc ${module} LOCATION)
    ADD_CUSTOM_COMMAND(
        TARGET vigranumpytest
        PRE_BUILD
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy_if_different ${loc} ../src/vigra/)
ENDFOREACH(module)

MACRO(DEPENDENCY_PATH variable target)
    GET_TARGET_PROPERTY(${variable} ${target} LOCATION)
    STRING(REGEX REPLACE "(/\\$\\([^\\)]*\\)/[^/]*|/[^/]*)$" "" ${variable} ${${variable}}) # get path prefix
    FILE(TO_NATIVE_PATH ${${variable}} ${variable})
ENDMACRO(DEPENDENCY_PATH)

DEPENDENCY_PATH(VIGRAIMPEX_PATH vigraimpex)
DEPENDENCY_PATH(VIGRANUMPYCMODULE_PATH vigranumpyore)
DEPENDENCY_PATH(VIGRANUMPYTEST_PATH vigranumpytest)

FILE(TO_CMAKE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../lib VIGRANUMPYSCRIPTS_PATH)
FILE(TO_NATIVE_PATH ${VIGRANUMPYSCRIPTS_PATH} VIGRANUMPYSCRIPTS_PATH)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/set_paths.py 
               ${CMAKE_CURRENT_BINARY_DIR}/set_paths.py 
               @ONLY)

IF(AUTOEXEC_TESTS)
    IF(NOT PYTHON_NOSETESTS_NOT_FOUND)
        IF (CMAKE_CFG_INTDIR STREQUAL ".")
            SET(TEST_ARGS -c "\"import" "nose\;" "nose.main()\"" .)
        ELSE()
            SET(TEST_ARGS -c "\"import" "nose\;" "nose.main()\"" . ${CMAKE_CFG_INTDIR})
        ENDIF()

        add_custom_command(
            TARGET vigranumpytest
            POST_BUILD
            COMMAND ${PYTHON_EXECUTABLE} ARGS ${TEST_ARGS}
            COMMENT "Running vigranumpy tests")
    ELSE()
        MESSAGE(STATUS "  vigranumpy tests will NOT be executed (nosetests missing)")
    ENDIF()
ENDIF()
