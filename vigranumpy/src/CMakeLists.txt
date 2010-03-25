ADD_SUBDIRECTORY(core)

IF(FFTW3_FOUND)
    ADD_SUBDIRECTORY(fourier)
ELSE()
    MESSAGE(STATUS "vigranumpy fourier module will not be built (fftw3 missing)")
ENDIF(FFTW3_FOUND)
