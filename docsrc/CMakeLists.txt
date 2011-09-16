IF(DOXYGEN_FOUND AND PYTHONINTERP_FOUND)
    execute_process(COMMAND ${DOXYGEN_EXECUTABLE} --version
                  OUTPUT_VARIABLE DOXYGEN_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

	CONFIGURE_FILE(
		${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
		${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		@ONLY)

    ADD_CUSTOM_TARGET (doc_cpp
		COMMAND ${CMAKE_COMMAND} -E make_directory ${DOCDIR}
		COMMAND ${CMAKE_COMMAND} -E make_directory ${DOCDIR}/vigra
		COMMAND ${DOXYGEN_EXECUTABLE}  
		COMMENT "Generating documentation with doxygen ${DOXYGEN_VERSION}")

	ADD_CUSTOM_COMMAND(
		TARGET doc_cpp
		POST_BUILD
		COMMAND ${PYTHON_EXECUTABLE}
			${CMAKE_CURRENT_SOURCE_DIR}/makeFunctionIndex.py
			${DOCDIR}/vigra
		COMMAND ${PYTHON_EXECUTABLE}
			${CMAKE_CURRENT_SOURCE_DIR}/post.py
			${DOCDIR}/vigra
			${vigra_version}
		COMMAND ${CMAKE_COMMAND} -E copy_directory
			${CMAKE_CURRENT_SOURCE_DIR}/documents
			${DOCDIR}/vigra/documents
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.txt
			${DOCDIR}/vigra/LICENSE.txt
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			${CMAKE_CURRENT_SOURCE_DIR}/../src/matlab/doc/tutorial_reference.pdf
			${DOCDIR}/vigra/documents/tutorial_reference.pdf
		COMMENT "Postprocessing html files")

ELSE()
	ADD_CUSTOM_TARGET (doc_cpp 
		COMMAND ${CMAKE_COMMAND} -E echo
			"Cannot generate C++ documentation "
			"(doxygen or Python not found)")

ENDIF(DOXYGEN_FOUND AND PYTHONINTERP_FOUND)

# create fallback index.html that simply refers to the online documentation
file(MAKE_DIRECTORY ${DOCDIR}/vigra)
file(GLOB INDEX_HTML_EXISTS ${DOCDIR}/vigra/index.html)
configure_file(vigra-icon.ico ${DOCDIR}/vigra/vigra-icon.ico COPYONLY)
if(NOT INDEX_HTML_EXISTS)
    configure_file(index_fallback.html ${DOCDIR}/vigra/index.html COPYONLY)
    configure_file(vigra.css ${DOCDIR}/vigra/vigra.css COPYONLY)
endif()

INSTALL(DIRECTORY ${DOCDIR}/
         DESTINATION ${DOCINSTALL}/)

ADD_DEPENDENCIES(doc doc_cpp)
