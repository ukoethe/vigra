ADD_SUBDIRECTORY(core)

#IF(FFTW3F_FOUND)
IF(0)
    ADD_SUBDIRECTORY(fourier)
ELSE()
    MESSAGE(STATUS "vigranumpy fourier module will not be built (fftw3f missing)")
ENDIF(FFTW3F_FOUND)
