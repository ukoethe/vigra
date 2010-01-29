IF(PYTHON_SPHINX)
	# vigranumpy documentation generation using sphinx
	SET(VIGRANUMPY_DOCDIR ${DOCDIR}/vigranumpy)
	SET(SPHINX_OPTS ${SPHINX_OPTS} CACHE STRING
		"options passed to sphinx during vigranumpy doc generation"
		FORCE)
	IF(NOT DEFINED SPHINX_PAPER)
		SET(SPHINX_PAPER a4)
	ENDIF()
	SET(SPHINX_PAPER ${SPHINX_PAPER} CACHE STRING
		"Papersize for sphinx documentation generation (a4, letter)"
		FORCE)
	SET(SPHINX_ALLOPTS
		-d ${VIGRANUMPY_DOCDIR}/doctrees
		-D latex_paper_size=${SPHINX_PAPER}
		${SPHINX_OPTS})
    
	ADD_CUSTOM_TARGET (doc_python
		COMMAND ${CMAKE_COMMAND} -E
			make_directory ${VIGRANUMPY_DOCDIR}
		COMMAND ${PYTHON_SPHINX} -b html
			${SPHINX_ALLOPTS}
			${CMAKE_CURRENT_SOURCE_DIR}
			${VIGRANUMPY_DOCDIR}/html
		COMMENT "Generating vigranumpy documentation")

ELSE(PYTHON_SPHINX)
	# no vigranumpy documentation if sphinx not available
	ADD_CUSTOM_TARGET (doc_python 
		${CMAKE_COMMAND} -E echo
		"Cannot generate Python documentation for vigranumpy. "
		"(sphinx-build not found)")

ENDIF(PYTHON_SPHINX)

ADD_DEPENDENCIES(doc doc_python)
