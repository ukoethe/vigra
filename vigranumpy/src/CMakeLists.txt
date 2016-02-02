VIGRA_CONFIGURE_THREADING()

ADD_SUBDIRECTORY(core)

IF(FFTW3F_FOUND)
    ADD_SUBDIRECTORY(fourier)
ELSE()
    MESSAGE(STATUS "vigranumpy fourier module will not be built (fftw missing)")
ENDIF()
