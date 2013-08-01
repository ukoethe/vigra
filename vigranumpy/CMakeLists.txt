IF(VIGRANUMPY_DEPENDENCIES_FOUND)
    ADD_CUSTOM_TARGET(vigranumpy)
    ADD_DEPENDENCIES(vigranumpy vigraimpex)

    SET(vigranumpy_tmp_dir "${CMAKE_BINARY_DIR}/vigranumpy/vigra")
    FILE(MAKE_DIRECTORY ${vigranumpy_tmp_dir})

    INCLUDE_DIRECTORIES(${VIGRANUMPY_INCLUDE_DIRS})

    # configure boost's autolink magic to use the right library name
    # (default on Windows is a mangled name like 'boost_python-vc110-mt-1_51.lib')
    if("${Boost_PYTHON_LIBRARY}" MATCHES "boost_python\\.lib")
        ADD_DEFINITIONS(-DBOOST_AUTO_LINK_NOMANGLE)
    endif()
        
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

    IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/private/CMakeLists.txt)
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
    
    set(VIGRANUMPY_CONFIG_INSTALL_PATH lib${LIB_SUFFIX}/vigranumpy)
    
    set(VIGRANUMPY_INCLUDE_RELATIVE "")
    foreach(PATH ${CMAKE_INSTALL_PREFIX}/include ${VIGRANUMPY_INCLUDE_DIRS})
        file(RELATIVE_PATH RELPATH ${CMAKE_INSTALL_PREFIX}/${VIGRANUMPY_CONFIG_INSTALL_PATH} ${PATH})
        set(VIGRANUMPY_INCLUDE_RELATIVE ${VIGRANUMPY_INCLUDE_RELATIVE} ${RELPATH})
    endforeach(PATH)
    
    set(VIGRANUMPY_LIBRARIES_RELATIVE "")
    foreach(PATH ${VIGRANUMPY_LIBRARIES})
        file(RELATIVE_PATH RELPATH ${CMAKE_INSTALL_PREFIX}/${VIGRANUMPY_CONFIG_INSTALL_PATH} ${PATH})
        set(VIGRANUMPY_LIBRARIES_RELATIVE ${VIGRANUMPY_LIBRARIES_RELATIVE} ${RELPATH})
    endforeach(PATH)
    
    CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/VigranumpyConfig.cmake.in
            ${PROJECT_BINARY_DIR}/lib/vigranumpy/CMake/VigranumpyConfig.cmake
            @ONLY IMMEDIATE)
    INSTALL(FILES ${PROJECT_BINARY_DIR}/lib/vigranumpy/CMake/VigranumpyConfig.cmake
            DESTINATION ${VIGRANUMPY_CONFIG_INSTALL_PATH})
ENDIF()


