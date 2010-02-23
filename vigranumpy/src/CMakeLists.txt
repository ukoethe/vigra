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
    ENDIF()
    
    if(DEFINED LIBRARIES)
        TARGET_LINK_LIBRARIES(${TARGET_NAME} ${LIBRARIES})
    endif()
    
    if(LIBRARY_NAME MATCHES "^vigranumpycore$")
        TARGET_LINK_LIBRARIES(${TARGET_NAME} ${VIGRANUMPY_LIBRARIES})
    else()
        TARGET_LINK_LIBRARIES(${TARGET_NAME} ${VIGRANUMPY_LIBRARIES} vigranumpy_core)
    endif()
    
    IF(PYTHON_PLATFORM MATCHES "^win32$")
        SET_TARGET_PROPERTIES(${TARGET_NAME} PROPERTIES OUTPUT_NAME "${LIBRARY_NAME}" 
                                                           PREFIX "" SUFFIX  ".pyd")
    ELSE()
        SET_TARGET_PROPERTIES(${TARGET_NAME} PROPERTIES OUTPUT_NAME "${LIBRARY_NAME}" 
                                                           PREFIX "")
    ENDIF()
    
    IF(PART_OF_VIGRANUMPY)
        IF(PYTHON_PLATFORM MATCHES "^win32$")
            INSTALL(TARGETS ${TARGET_NAME} RUNTIME DESTINATION ${VIGRANUMPY_INSTALL_DIR}/vigra)
        ELSE()
            INSTALL(TARGETS ${TARGET_NAME}
                    LIBRARY DESTINATION ${VIGRANUMPY_INSTALL_DIR}/vigra)
        ENDIF()
        
        # create a temporary vigranumpy installation in ${vigranumpy_tmp_dir}
        # (required for testing and documentation generation)
        GET_TARGET_PROPERTY(loc ${TARGET_NAME} LOCATION)
        ADD_CUSTOM_COMMAND(
            TARGET ${TARGET_NAME}
            POST_BUILD
            COMMAND ${CMAKE_COMMAND}
            ARGS -E copy_if_different ${loc} ${vigranumpy_tmp_dir}/)
    ENDIF()
ENDFUNCTION(VIGRA_ADD_NUMPY_MODULE)

###############################################################################

INCLUDE_DIRECTORIES(${VIGRANUMPY_INCLUDE_DIRS})

IF(HDF5_FOUND)
    ADD_DEFINITIONS(${HDF5_CPPFLAGS})
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
ENDIF(HDF5_FOUND)

IF(MSVC OR MINGW)
    ADD_DEFINITIONS(-DVIGRA_DLL)
ENDIF ()

ADD_SUBDIRECTORY(core)

IF(FFTW3_FOUND)
    ADD_SUBDIRECTORY(fourier)
ELSE()
    MESSAGE(STATUS "vigranumpy fourier module will not be built (fftw3 missing)")
ENDIF(FFTW3_FOUND)

ADD_SUBDIRECTORY(experiments EXCLUDE_FROM_ALL)
