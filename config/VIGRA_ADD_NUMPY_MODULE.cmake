###################################################################
#
# VIGRA_ADD_NUMPY_MODULE: setup a module dependening on vigranumpy
#
# VIGRA_ADD_NUMPY_MODULE(modulename [SOURCES] source1.cpp, source2.cpp ...
#                                   [LIBRARIES dependency1 dependency2 ...]
#                                   [VIGANUMPY])
#
#        'modulename' is the module name to be used within Python (e.g. 'import modulename').
#        Unless 'VIGRANUMPY' is specified (see below), it is also the cmake target name.
#
#        SOURCE are the C++ sources of the module, LIBRARIES the necessary libraries.
#        Dependency syntax must conform to the requirements of the cmake command
#        TARGET_LINK_LIBRARIES. Modules are automatically linked against vigranumpycore
#        and its dependencies (libpython, boost_python), so it is not necessary to state
#        this dependency explicitly.
#
#        If VIGRANUMPY is given, the module is considered part of 'vigranumpy' and will
#        be compiled and installed along with the other vigranumpy modules (otherwise,
#        no installation target will be defined). The cmake target name becomes
#        'vigranumpy_modulename' in order to get useful alphabetic sorting of
#        targets in project files.
FUNCTION(VIGRA_ADD_NUMPY_MODULE target)

    # parse the args
    set(v SOURCES)
    set(PART_OF_VIGRANUMPY 0)
    foreach(i ${ARGN})
        if(${i} MATCHES "^SOURCES$")
            set(v SOURCES)
        elseif(${i} MATCHES "^LIBRARIES$")
            set(v LIBRARIES)
        elseif(${i} MATCHES "^VIGRANUMPY$")
            set(PART_OF_VIGRANUMPY 1)
        else()
            set(${v} ${${v}} ${i})
        endif()
    endforeach(i)

    IF(PART_OF_VIGRANUMPY)
        set(TARGET_NAME vigranumpy_${target})
        if(target MATCHES "^core$")
            set(LIBRARY_NAME vigranumpycore)
        else()
            set(LIBRARY_NAME ${target})
        endif()
    ELSE()
        set(TARGET_NAME ${target})
        set(LIBRARY_NAME ${target})
    ENDIF()

    ADD_LIBRARY(${TARGET_NAME} SHARED ${SOURCES})

    IF(PART_OF_VIGRANUMPY)
        ADD_DEPENDENCIES(vigranumpy ${TARGET_NAME})

        # Store dependencies as a custom target property, so that we can
        # later query them.
        # TODO: Does cmake provide a standard way to query the dependencies?
        GET_TARGET_PROPERTY(VIGRANUMPY_DEPENDS vigranumpy VIGRA_DEPENDS)
        IF(NOT VIGRANUMPY_DEPENDS)
            set(VIGRANUMPY_DEPENDS "")
        ENDIF()
        SET_TARGET_PROPERTIES(vigranumpy PROPERTIES VIGRA_DEPENDS "${VIGRANUMPY_DEPENDS} ${TARGET_NAME}")
    ENDIF()

    ######################################################################
    #
    #      find Python EXT SUFFIX for compatibility with PyPy and CPython
    #
    ######################################################################
    execute_process(COMMAND ${Python_EXECUTABLE} -c
                     "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
                     OUTPUT_VARIABLE EXT_SUFFIX OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(DEFINED LIBRARIES)
        TARGET_LINK_LIBRARIES(${TARGET_NAME} ${LIBRARIES})
    endif()

    TARGET_LINK_LIBRARIES(${TARGET_NAME} ${VIGRANUMPY_LIBRARIES})

    IF(PYTHON_PLATFORM MATCHES "^windows$")
        SET_TARGET_PROPERTIES(${TARGET_NAME} PROPERTIES OUTPUT_NAME "${LIBRARY_NAME}"
                                                           PREFIX "" SUFFIX  ${EXT_SUFFIX})
    ELSEIF(MACOSX)
        IF(${CMAKE_MAJOR_VERSION} LESS 3)
            SET_TARGET_PROPERTIES(${TARGET_NAME} PROPERTIES OUTPUT_NAME "${LIBRARY_NAME}" PREFIX ""
                              SUFFIX ${EXT_SUFFIX} INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${VIGRANUMPY_INSTALL_DIR}/vigra")
        ELSE()
            SET_TARGET_PROPERTIES(${TARGET_NAME} PROPERTIES OUTPUT_NAME "${LIBRARY_NAME}" PREFIX ""    SUFFIX ${EXT_SUFFIX} )
        ENDIF()
    ELSE()
        SET_TARGET_PROPERTIES(${TARGET_NAME} PROPERTIES OUTPUT_NAME "${LIBRARY_NAME}"
                                                           PREFIX "" SUFFIX ${EXT_SUFFIX})
    ENDIF()

    IF(PART_OF_VIGRANUMPY)
        IF(PYTHON_PLATFORM MATCHES "^windows$")
            INSTALL(TARGETS ${TARGET_NAME} RUNTIME DESTINATION ${VIGRANUMPY_INSTALL_DIR}/vigra)
        ELSE()
            INSTALL(TARGETS ${TARGET_NAME}
                    LIBRARY DESTINATION ${VIGRANUMPY_INSTALL_DIR}/vigra)
        ENDIF()

        # create a temporary vigranumpy installation in ${vigranumpy_tmp_dir}
        # (required for testing and documentation generation)
        ADD_CUSTOM_COMMAND(
            TARGET ${TARGET_NAME}
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} ARGS -E copy_if_different $<TARGET_FILE:${TARGET_NAME}> ${vigranumpy_tmp_dir}/
            COMMENT "Copying target ${TARGET_NAME} to temporary module directory")
    ENDIF()
ENDFUNCTION(VIGRA_ADD_NUMPY_MODULE)
