SET(PYSOURCES
    __init__.py
    imagewindow.py
    quickdialog.py
    )

INSTALL(FILES ${PYSOURCES} DESTINATION ${VIGRANUMPY_INSTALL_DIR}/vigra/pyqt)

FOREACH(lib_file ${PYSOURCES})
    ADD_CUSTOM_COMMAND(
        TARGET vigranumpy_lib
        POST_BUILD
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${lib_file} ${vigranumpy_tmp_dir}/pyqt/${lib_file}
        COMMENT "Copying Python sources to temporary module directory")
ENDFOREACH(lib_file)
