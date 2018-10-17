INCLUDE_DIRECTORIES(${SUPPRESS_WARNINGS} ${VIGRANUMPY_INCLUDE_DIRS})

SET(TEST_SCRIPTS
    test_arraytypes.py
    test_rf.py
    test_impex.py
    test1.py
    test2.py
    test3.py
    test4.py
    test_color.py
    test_segmentation.py
    test_multidef.py
    )

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
IF(AUTOBUILD_TESTS)
  ADD_LIBRARY(vigranumpytest SHARED
               vigranumpytest.cxx ${TEST_SUCCESS_FILE} ${TEST_SCRIPTS})
ELSE()
  ADD_LIBRARY(vigranumpytest SHARED EXCLUDE_FROM_ALL
               vigranumpytest.cxx ${TEST_SUCCESS_FILE} ${TEST_SCRIPTS})
ENDIF()

IF(PYTHON_PLATFORM MATCHES "^windows$")
    SET_TARGET_PROPERTIES(vigranumpytest PROPERTIES PREFIX "" SUFFIX  ".pyd")
ELSEIF(MACOSX)
    SET_TARGET_PROPERTIES(vigranumpytest PROPERTIES PREFIX "" SUFFIX ".so")
ELSE()
    SET_TARGET_PROPERTIES(vigranumpytest PROPERTIES PREFIX "")
ENDIF()

TARGET_LINK_LIBRARIES(vigranumpytest ${VIGRANUMPY_LIBRARIES})
ADD_DEPENDENCIES(check_python vigranumpytest)
ADD_DEPENDENCIES(ctest vigranumpytest)
ADD_DEPENDENCIES(vigranumpytest vigranumpy)

# setup the file 'set_path.py' to ensure that tests are run against the binaries
# which are currently being built, not against (possibly outdated) installed modules
STRING(REGEX REPLACE "/vigra$" "" VIGRANUMPY_TMP_PATH ${vigranumpy_tmp_dir})
FILE(TO_NATIVE_PATH "${VIGRANUMPY_TMP_PATH}" VIGRANUMPY_TMP_PATH)

if(CMAKE_MAJOR_VERSION LESS 3)
    DEPENDENCY_PATH(VIGRAIMPEX_PATH vigraimpex)
    DEPENDENCY_PATH(VIGRANUMPYTEST_PATH vigranumpytest)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/set_paths.py.cmake2.in
                   ${CMAKE_CURRENT_BINARY_DIR}/set_paths.py
                   @ONLY)
else()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/set_paths.py.in
                   ${CMAKE_CURRENT_BINARY_DIR}/set_paths.py.in
                   @ONLY)

    # two-stage file configuration is necessary because certain target
    # properties are only known at generation time (policy CMP0026)
    if(MSVC)
        file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/set_paths.py
                      INPUT  ${CMAKE_CURRENT_BINARY_DIR}/set_paths.py.in
                      CONDITION $<CONFIG:Release>)
    else()
        file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/set_paths.py
                      INPUT  ${CMAKE_CURRENT_BINARY_DIR}/set_paths.py.in)
    endif()
endif()

IF(NOT PYTHON_NOSETESTS_NOT_FOUND)
    # copy the individual test scripts
    FOREACH(test_script ${TEST_SCRIPTS})
        ADD_CUSTOM_COMMAND(
            TARGET vigranumpytest
            POST_BUILD
            COMMAND ${CMAKE_COMMAND}
            ARGS -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/${test_script} ${CMAKE_CURRENT_BINARY_DIR}/${test_script}
            COMMENT "Running vigranumpy tests")
    ENDFOREACH(test_script)

    VIGRA_NATIVE_PATH(PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE})

    # Windows: Set the DLL path to the binaries currently being built
    #          (so that binaries which are already installed are ignored)
    #          This is not necessary on Linux, because the correct path
    #          is hard-coded into the module by means of 'rpath'
    set(EXTRA_PATH "")
    IF(MSVC)
        SET(PATHSEP ";")
    ELSE()
        SET(PATHSEP ":")
    ENDIF()
    FOREACH(lib ${LIBRARIES})
        if(TARGET ${lib})
            SET(EXTRA_PATH  "$<TARGET_FILE_DIR:${lib}>${PATHSEP}${EXTRA_PATH}")
        endif()
    ENDFOREACH(lib)

    # Setup the test execution script
    # two-stage file configuration is necessary because certain target
    # properties are only known at generation time (policy CMP0026)
    VIGRA_NATIVE_PATH(VIGRA_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
    IF(MSVC OR MINGW)
        get_filename_component(PYTHON_PATH ${PYTHON_EXECUTABLE} PATH)
        get_filename_component(PYTHON_EXE  ${PYTHON_EXECUTABLE} NAME)
        SET(EXTRA_PATH  "${PYTHON_PATH}${PATHSEP}${EXTRA_PATH}")
        SET(VIGRA_TEST_EXECUTABLE "${PYTHON_EXE} -c \"import nose; nose.main()\" . %CONFIGURATION%")
        SET(VIGRA_TEST_SCRIPT     "${CMAKE_CURRENT_BINARY_DIR}/run_vigranumpytest.bat")
        CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/config/run_test.bat.in
                       ${VIGRA_TEST_SCRIPT}.in
                       @ONLY)
        file(GENERATE OUTPUT "${VIGRA_TEST_SCRIPT}"
                      INPUT  "${VIGRA_TEST_SCRIPT}.in")
    ELSE()
        SET(VIGRA_TEST_EXECUTABLE "${PYTHON_EXECUTABLE} -c \"import nose; nose.main()\" .")
        SET(VIGRA_TEST_SCRIPT     "${CMAKE_CURRENT_BINARY_DIR}/run_vigranumpytest.sh")
        CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/config/run_test.sh.in
                       ${VIGRA_TEST_SCRIPT}.in
                       @ONLY)
        EXECUTE_PROCESS(COMMAND chmod u+x ${VIGRA_TEST_SCRIPT}.in OUTPUT_QUIET ERROR_QUIET)
        file(GENERATE OUTPUT "${VIGRA_TEST_SCRIPT}"
                      INPUT  "${VIGRA_TEST_SCRIPT}.in")
    ENDIF()

    # register the test execution command
    # VIGRA_CONFIGURATION is 'release' or 'debug' on Windows, nothing on Linux
    # variable CMAKE_CFG_INTDIR contains a dot '.' for a Windows nmake build, or
    #                           '$(OutDir)' for a VisualStudio build (OutDir will be
    #                           set by VS at build time)
    SET(VIGRA_CONFIGURATION ${CMAKE_CFG_INTDIR})

    IF(AUTOEXEC_TESTS)
        add_custom_command(
            TARGET vigranumpytest
            POST_BUILD
            COMMAND ${VIGRA_TEST_SCRIPT} ARGS ${VIGRA_CONFIGURATION}
            COMMENT "Running vigranumpy tests")
    ENDIF()

    ADD_TEST(vigranumpytest ${VIGRA_TEST_SCRIPT} ${VIGRA_CONFIGURATION})

    # if we configure for Visual Studio, setup the debug command
    if(MSVC AND NOT (MSVC_VERSION VERSION_LESS "11"))
        # FIXME: this may not portable between VC versions (works for 11.0 - 14.0)
        set(VIGRA_TEST_EXE "${PYTHON_EXECUTABLE}")
        set(VIGRA_TEST_ARGS "-c \"import nose; nose.main()\" . \$(Configuration)")
        configure_file(${CMAKE_SOURCE_DIR}/config/testdebug.vcxproj.user.in
                       ${CMAKE_CURRENT_BINARY_DIR}/vigranumpytest.vcxproj.user
                       @ONLY)
    endif()
ELSE()
        MESSAGE(STATUS "  vigranumpy tests will NOT be executed (nosetests missing)")
ENDIF()
