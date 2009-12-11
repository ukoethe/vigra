IF(PYTHON_SPHINX)
    SET(VIGRANUMPY_DOCDIR ${vigra_SOURCE_DIR}/vigranumpy/doc)
    
    FILE(MAKE_DIRECTORY ${VIGRANUMPY_DOCDIR})

    if(MSVC OR MINGW)
        add_custom_target (doc_python ${CMAKE_CURRENT_SOURCE_DIR}/make.bat html 
                               WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
                               COMMENT "Generating vigranumpy documentation")
    else()
        add_custom_target (doc_python make html
                               WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
                               COMMENT "Generating vigranumpy documentation")
    endif()
    ADD_DEPENDENCIES(doc doc_python)

    INSTALL(DIRECTORY ${VIGRANUMPY_DOCDIR} 
            DESTINATION ${VIGRANUMPY_INSTALL_DIR}/vigra)

ELSE()
    add_custom_target (doc_python 
       echo "Cannot generate Python documentation (sphinx-build not found)")

    ADD_DEPENDENCIES(doc doc_python)

ENDIF()
