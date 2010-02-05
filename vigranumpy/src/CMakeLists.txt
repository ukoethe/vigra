ADD_SUBDIRECTORY(core vigra)

IF(FFTW3_FOUND)
    ADD_SUBDIRECTORY(fourier)
ELSE()
    MESSAGE(STATUS "vigranumpy fourier module will not be built (fftw3 missing)")
ENDIF(FFTW3_FOUND)

ADD_CUSTOM_TARGET(build_dir_vigranumpy)
ADD_DEPENDENCIES(build_dir_vigranumpy vigranumpy)

SET(tmp_vigranumpy_dir ${CMAKE_BINARY_DIR}/vigranumpy/src/vigra)

#copy rest file into build-vigra module
SET(LIB_FILES __init__.py arraytypes.py ufunc.py)
FOREACH(lib_file ${LIB_FILES})
    ADD_CUSTOM_COMMAND(
        TARGET build_dir_vigranumpy
        PRE_BUILD
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy_if_different
        ${PROJECT_SOURCE_DIR}/vigranumpy/lib/${lib_file} ${tmp_vigranumpy_dir}/${lib_file}
        )
ENDFOREACH(lib_file)

#copy vigranumpy files into build-vigra module
SET(MODULE_TARGETS vigranumpycore convolution segmentation edgedetection
    morphology classification tensor impex filters noise)

FOREACH(module ${MODULE_TARGETS})
    GET_TARGET_PROPERTY(loc ${module} LOCATION)
    ADD_CUSTOM_COMMAND(
        TARGET build_dir_vigranumpy
        PRE_BUILD
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy_if_different ${loc} ${tmp_vigranumpy_dir}/)
ENDFOREACH(module)
