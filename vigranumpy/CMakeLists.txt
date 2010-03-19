IF(VIGRANUMPY_DEPENDENCIES_FOUND)
    ADD_CUSTOM_TARGET(vigranumpy)
    ADD_DEPENDENCIES(vigranumpy vigraimpex)

    SET(vigranumpy_tmp_dir "${CMAKE_BINARY_DIR}/vigranumpy/vigra")
    FILE(MAKE_DIRECTORY ${vigranumpy_tmp_dir})

    ADD_SUBDIRECTORY(src)
    ADD_SUBDIRECTORY(lib)
    ADD_SUBDIRECTORY(test)
    ADD_SUBDIRECTORY(docsrc)
ENDIF()
