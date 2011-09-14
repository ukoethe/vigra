FILE(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/__version__.py
     "version = '${vigra_version}'
")
SET(PYSOURCES
    __init__.py
    __version__.py
    tagged_array.py
    arraytypes.py
    ufunc.py
    )

INSTALL(FILES ${PYSOURCES} DESTINATION ${VIGRANUMPY_INSTALL_DIR}/vigra)

ADD_CUSTOM_TARGET(vigranumpy_lib)
ADD_DEPENDENCIES(vigranumpy vigranumpy_lib)

FOREACH(lib_file ${PYSOURCES})
    ADD_CUSTOM_COMMAND(
        TARGET vigranumpy_lib
        POST_BUILD
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${lib_file} ${vigranumpy_tmp_dir}/${lib_file}
        COMMENT "Copying Python sources to temporary module directory")
ENDFOREACH(lib_file)
ADD_CUSTOM_COMMAND(
    TARGET vigranumpy_lib
    POST_BUILD
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pyqt ${vigranumpy_tmp_dir}/pyqt
    COMMENT "Copying Python sources to temporary module directory")

ADD_SUBDIRECTORY(pyqt)