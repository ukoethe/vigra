IF(DOXYGEN_FOUND AND PYTHONINTERP_FOUND)
	CONFIGURE_FILE(
		${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
		${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		@ONLY)

    ADD_CUSTOM_TARGET (doc_cpp
		COMMAND ${CMAKE_COMMAND} -E make_directory ${DOCDIR}
		COMMAND ${CMAKE_COMMAND} -E make_directory ${DOCDIR}/vigra
		COMMAND ${DOXYGEN_EXECUTABLE}  
		COMMENT "Generating doxygen documentation")

	ADD_CUSTOM_COMMAND(
		TARGET doc_cpp
		POST_BUILD
		COMMAND ${PYTHON_EXECUTABLE}
			${CMAKE_CURRENT_SOURCE_DIR}/makeFunctionIndex.py
			${DOCDIR}/vigra
		COMMAND ${PYTHON_EXECUTABLE}
			${CMAKE_CURRENT_SOURCE_DIR}/post-doxy-1.5.6.py
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

	INSTALL(DIRECTORY ${DOCDIR}/
		DESTINATION ${DOCINSTALL}/)

ELSE()
	ADD_CUSTOM_TARGET (doc_cpp 
		COMMAND ${CMAKE_COMMAND} -E echo
			"Cannot generate C++ documentation "
			"(doxygen or Python not found)")

ENDIF(DOXYGEN_FOUND AND PYTHONINTERP_FOUND)

ADD_DEPENDENCIES(doc doc_cpp)
