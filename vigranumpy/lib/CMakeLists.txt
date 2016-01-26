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

SET(TIMESTAMP "${CMAKE_CURRENT_BINARY_DIR}/TIMESTAMP.cxx")
FILE(GLOB TIMESTAMP_FOUND ${TIMESTAMP})
IF(NOT TIMESTAMP_FOUND)
    FILE(WRITE ${TIMESTAMP}
     "// auto-generated dummy file to ensure dependency on Python files. 
// ADD_CUSTOM_TARGET(vigranumpy_lib SOURCES 'dollarsign'{PYSOURCES}) does not seem to be sufficient.
")
ENDIF()

add_custom_command(
    OUTPUT ${TIMESTAMP}
    DEPENDS ${PYSOURCES}
    COMMAND ${CMAKE_COMMAND}
    ARGS -E touch ${TIMESTAMP})

ADD_CUSTOM_TARGET(vigranumpy_lib SOURCES ${TIMESTAMP} ${PYSOURCES})
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