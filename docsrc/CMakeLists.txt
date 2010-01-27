IF(DOXYGEN_FOUND AND PYTHONINTERP_FOUND)

    FILE(MAKE_DIRECTORY ${DOCDIR})

    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

    add_custom_target (doc_cpp ${DOXYGEN_EXECUTABLE}  
                               COMMENT "Generating documentation")

    ADD_DEPENDENCIES(doc doc_cpp)

    add_custom_command(
            TARGET doc_cpp
            POST_BUILD
            COMMAND ${PYTHON_EXECUTABLE} ARGS makeFunctionIndex.py ${DOCDIR}
            COMMAND ${PYTHON_EXECUTABLE} ARGS post-doxy-1.5.6.py ${DOCDIR} ${vigra_version}
            COMMAND ${CMAKE_COMMAND} 
               ARGS -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/documents ${DOCDIR}/documents
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
            COMMENT "Postprocessing html files")

    INSTALL(DIRECTORY ${DOCDIR}/ DESTINATION ${DOCINSTALL})

ELSE()
    add_custom_target (doc_cpp 
       echo "Cannot generate C++ documentation (doxygen or Python not found)")

    ADD_DEPENDENCIES(doc doc_cpp)

ENDIF()
