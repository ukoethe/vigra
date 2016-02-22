IF(PYTHON_SPHINX)
    # vigranumpy documentation generation using sphinx
    IF(DEFINED VIGRANUMPY_DOCDIR)
        MESSAGE(SEND_ERROR "VIGRANUMPY_DOCDIR already defined! Please delete from Cache.")
    ENDIF()
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
        -c ${CMAKE_CURRENT_BINARY_DIR}
        ${SPHINX_OPTS})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/_static)
    configure_file(_static/vigra-icon.ico ${CMAKE_CURRENT_BINARY_DIR}/_static/vigra-icon.ico COPYONLY)

    STRING(REGEX REPLACE "/vigra$" "" VIGRANUMPY_TMP_PATH ${vigranumpy_tmp_dir})

    if(CMAKE_MAJOR_VERSION LESS 3)
        DEPENDENCY_PATH(VIGRAIMPEX_PATH vigraimpex)
        CONFIGURE_FILE(
            ${CMAKE_CURRENT_SOURCE_DIR}/conf.py.cmake2.in
            ${CMAKE_CURRENT_BINARY_DIR}/conf.py
            @ONLY)
    else()
        CONFIGURE_FILE(
            ${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in
            ${CMAKE_CURRENT_BINARY_DIR}/conf.py.in
            @ONLY)

        # two-stage file configuration is necessary because certain target
        # properties are only known at generation time (policy CMP0026)
        if(MSVC)
            file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/conf.py
                          INPUT  ${CMAKE_CURRENT_BINARY_DIR}/conf.py.in
                          CONDITION $<CONFIG:Release>)
        else()
            file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/conf.py
                          INPUT  ${CMAKE_CURRENT_BINARY_DIR}/conf.py.in)
        endif()
    endif()

    # VIGRA_CONFIGURATION is 'release' or 'debug' on Windows, nothing on Linux
    # variable CMAKE_CFG_INTDIR contains a dot '.' for a Windows nmake build, or
    #                           '$(OutDir)' for a VisualStudio build (OutDir will be
    #                           set by VS at build time)
    IF(CMAKE_CFG_INTDIR AND NOT CMAKE_CFG_INTDIR STREQUAL ".")
        SET(VIGRA_CONFIGURATION -D build_configuration=${CMAKE_CFG_INTDIR})
    ELSE()
        SET(VIGRA_CONFIGURATION)
    ENDIF()

    # just rerun sphinx (useful to debug the vigranumpy docu)
    ADD_CUSTOM_TARGET (doc_sphinx
        COMMAND ${CMAKE_COMMAND} -E make_directory
            ${VIGRANUMPY_DOCDIR}
        COMMAND ${PYTHON_SPHINX} -b html
            ${SPHINX_ALLOPTS}
            ${VIGRA_CONFIGURATION}
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${VIGRANUMPY_DOCDIR}
        COMMENT "Generating vigranumpy documentation")

    # create vigranumpy documentation after updating C++ docu and vigranumpy modules
    ADD_CUSTOM_TARGET(doc_python
        COMMAND ${CMAKE_COMMAND} -E make_directory
            ${VIGRANUMPY_DOCDIR}
        COMMAND ${PYTHON_SPHINX} -b html
            ${SPHINX_ALLOPTS}
            ${VIGRA_CONFIGURATION}
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${VIGRANUMPY_DOCDIR}
        COMMENT "Generating vigranumpy documentation")
    ADD_DEPENDENCIES(doc_python vigranumpy)
    ADD_DEPENDENCIES(doc_python doc_cpp)

ELSE(PYTHON_SPHINX)
    # no vigranumpy documentation if sphinx not available
    ADD_CUSTOM_TARGET (doc_python
        ${CMAKE_COMMAND} -E echo
        "Cannot generate Python documentation for vigranumpy. "
        "(sphinx-build not found)"
        VERBATIM)

ENDIF(PYTHON_SPHINX)

ADD_DEPENDENCIES(doc doc_python)
