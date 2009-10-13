# Enhanced version of ADD_TEST
#
# Usage of this module:
#
#     VIGRA_ADD_TEST(target [test_data1 test_data2 ...])
#
# This macro calls ADD_TEST(${target} ${target}) and additionally
# * copies the given test data files into the target directory (instead of interpreting
#   them as arguments to the test program, as ADD_TEST() would do)
# * installs a post-build event that runs the test automatically after linking (Visual Studio only)
#
INCLUDE(TestOrDelete RESULT_VARIABLE TEST_OR_DELETE)

IF (CMAKE_GENERATOR MATCHES "Visual Studio")
   STRING(REGEX REPLACE "\\.cmake$" ".bat" TEST_OR_DELETE ${TEST_OR_DELETE})
ELSE ()
   STRING(REGEX REPLACE "\\.cmake$" ".sh" TEST_OR_DELETE ${TEST_OR_DELETE})
ENDIF ()

MACRO(VIGRA_ADD_TEST target)
    ADD_DEPENDENCIES(test ${target})
    GET_TARGET_PROPERTY(${target}_executable ${target} LOCATION)
    file(TO_NATIVE_PATH ${${target}_executable} ${target}_executable)
    add_custom_command(
        TARGET ${target}
        POST_BUILD
        COMMAND ${TEST_OR_DELETE} ARGS ${${target}_executable}
        COMMENT "Running tests")
    FOREACH(test_data ${ARGN})
      add_custom_command(
        TARGET ${target}
        PRE_BUILD
        COMMAND ${CMAKE_COMMAND} ARGS -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${test_data} ${CMAKE_CURRENT_BINARY_DIR}/${test_data}
        COMMENT "Copying test data")
    ENDFOREACH(test_data)
ENDMACRO(VIGRA_ADD_TEST)

