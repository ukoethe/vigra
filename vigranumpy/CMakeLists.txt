IF(VIGRANUMPY_DEPENDENCIES_FOUND)
    ADD_CUSTOM_TARGET(vigranumpy)
    ADD_DEPENDENCIES(vigranumpy vigraimpex)

    SET(vigranumpy_tmp_dir "${CMAKE_BINARY_DIR}/vigranumpy/vigra")
    FILE(MAKE_DIRECTORY ${vigranumpy_tmp_dir})
    FILE(WRITE ${vigranumpy_tmp_dir}/__version__.py
         "version = '${vigra_version}'
")

    INCLUDE_DIRECTORIES(${VIGRANUMPY_INCLUDE_DIRS})

    IF(HDF5_FOUND)
        ADD_DEFINITIONS(-DHasHDF5 ${HDF5_CPPFLAGS})
        INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
        SET(VIGRANUMPY_IMPEX_LIBRARIES vigraimpex ${HDF5_LIBRARIES})
    ELSE()
        SET(VIGRANUMPY_IMPEX_LIBRARIES vigraimpex)
    ENDIF(HDF5_FOUND)

    IF(MSVC OR MINGW)
        ADD_DEFINITIONS(-DVIGRA_DLL)
    ENDIF ()
    
    include(VIGRA_ADD_NUMPY_MODULE)

    ADD_SUBDIRECTORY(src)
    ADD_SUBDIRECTORY(lib)
    ADD_SUBDIRECTORY(test)
    ADD_SUBDIRECTORY(docsrc)

    IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/private)
        ADD_SUBDIRECTORY(private EXCLUDE_FROM_ALL)
    ENDIF()
    
    IF(PYTHON_PLATFORM MATCHES "^windows$")
        FILE(MAKE_DIRECTORY ${vigranumpy_tmp_dir}/dlls)    
        GET_TARGET_PROPERTY(VIGRAIMPEX_LOCATION vigraimpex LOCATION)
        file(RELATIVE_PATH VIGRAIMPEX_DOCDIR ${vigranumpy_tmp_dir} ${DOCDIR})
        
        CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in
                        ${CMAKE_CURRENT_BINARY_DIR}/setup.py
                        @ONLY)

        add_custom_target(PACKAGE_VIGRANUMPY 
                           COMMAND ${CMAKE_COMMAND} -E copy_directory ${DOCDIR} ${vigranumpy_tmp_dir}/doc
                           COMMAND python setup.py bdist_wininst
                           DEPENDS check
                           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                           COMMENT "Creating vigranumpy setup in ${CMAKE_CURRENT_BINARY_DIR}/dist")        
    ENDIF()
ENDIF()


