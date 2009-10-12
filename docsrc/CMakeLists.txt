FIND_PACKAGE(Doxygen)
FIND_PACKAGE(PythonInterp)

SET(DOCDIR ${vigra_SOURCE_DIR}/doc/vigra)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile @ONLY)

add_custom_target (doc doxygen  
                       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
                       COMMENT "Generating documentation")

IF(WIN32)
    file(TO_NATIVE_PATH ${DOCDIR} NATIVE_DOCDIR)
    SET(DELETE_COMMAND del /F /S /Q)
ELSE()
    SET(NATIVE_DOCDIR ${DOCDIR})
    SET(DELETE_COMMAND rm -rf)
ENDIF()

add_custom_command(
        TARGET doc
        PRE_BUILD
        COMMAND ${DELETE_COMMAND} ARGS ${NATIVE_DOCDIR}
        COMMENT "Removing old doc files")

add_custom_command(
        TARGET doc
        POST_BUILD
        COMMAND ${PYTHON_EXECUTABLE} ARGS makeFunctionIndex.py ${DOCDIR}
        COMMAND ${PYTHON_EXECUTABLE} ARGS post-doxy-1.5.6.py ${DOCDIR} ${vigra_version}
        COMMAND ${CMAKE_COMMAND} ARGS -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/documents ${DOCDIR}/documents
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        COMMENT "Postprocessing html files")
