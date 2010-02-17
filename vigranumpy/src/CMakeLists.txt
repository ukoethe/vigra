FUNCTION(VIGRA_ADD_NUMPY_MODULE target)
    set(TARGET_NAME vigranumpy_${target})
    if(target MATCHES "^core$")
        set(LIBRARY_NAME vigranumpycore)
    else()
        set(LIBRARY_NAME ${target})
    endif()
    
    # parse the args
    set(v SOURCES)
    foreach(i ${ARGN})
        if(${i} MATCHES "^SOURCES$")
            set(v SOURCES)
        elseif(${i} MATCHES "^LIBRARIES$")
            set(v LIBRARIES)
        else()
            set(${v} ${${v}} ${i})
        endif()
    endforeach(i)
    
    ADD_LIBRARY(${TARGET_NAME} SHARED ${SOURCES})    
    ADD_DEPENDENCIES(vigranumpy ${TARGET_NAME})
    if(DEFINED LIBRARIES)
        TARGET_LINK_LIBRARIES(${TARGET_NAME} ${LIBRARIES})
    else()
        TARGET_LINK_LIBRARIES(${TARGET_NAME} ${VIGRANUMPY_LIBRARIES} vigranumpy_core)
    endif()
    
    IF(PYTHON_PLATFORM MATCHES "^win32$")
        SET_TARGET_PROPERTIES(${TARGET_NAME} PROPERTIES OUTPUT_NAME "${LIBRARY_NAME}" 
                                                           PREFIX "" SUFFIX  ".pyd")
        INSTALL(TARGETS ${TARGET_NAME} RUNTIME DESTINATION ${VIGRANUMPY_INSTALL_DIR}/vigra)
    ELSE()
        SET_TARGET_PROPERTIES(${TARGET_NAME} PROPERTIES OUTPUT_NAME "${LIBRARY_NAME}" 
                                                           PREFIX "")
        INSTALL(TARGETS ${TARGET_NAME}
                 LIBRARY DESTINATION ${VIGRANUMPY_INSTALL_DIR}/vigra)
    ENDIF()
    
    GET_TARGET_PROPERTY(loc ${TARGET_NAME} LOCATION)
    ADD_CUSTOM_COMMAND(
        TARGET ${TARGET_NAME}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy_if_different ${loc} ${vigranumpy_tmp_dir}/)
ENDFUNCTION(VIGRA_ADD_NUMPY_MODULE)

###############################################################################

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
