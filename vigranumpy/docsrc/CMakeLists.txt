IF(PYTHON_SPHINX)
    SET(VIGRANUMPY_DOCDIR ${DOCDIR}/vigranumpy)
	SET(SPHINX_OPTS)
	SET(SPHINX_PAPER a4)
	SET(SPHINX_PAPEROPT_a4 -D latex_paper_size=a4)
	SET(SPHINX_PAPEROPT_letter -D latex_paper_size=letter)
	SET(SPHINX_ALLOPTS -d ${VIGRANUMPY_DOCDIR}/doctrees ${SPHINX_PAPEROPT_${SPHINX_PAPER}} ${SPHINX_OPTS} ${CMAKE_CURRENT_SOURCE_DIR})
    
    FILE(MAKE_DIRECTORY ${VIGRANUMPY_DOCDIR})

    ADD_CUSTOM_TARGET (doc_python ${PYTHON_SPHINX} -b html ${SPHINX_ALLOPTS} ${VIGRANUMPY_DOCDIR}/html
                       COMMENT "Generating vigranumpy documentation")
    ADD_DEPENDENCIES(doc doc_python)

    INSTALL(DIRECTORY ${VIGRANUMPY_DOCDIR}/
            DESTINATION ${VIGRANUMPY_INSTALL_DIR}/vigra/doc)

ELSE()
    add_custom_target (doc_python 
       echo "Cannot generate Python documentation (sphinx-build not found)")

    ADD_DEPENDENCIES(doc doc_python)

ENDIF()
