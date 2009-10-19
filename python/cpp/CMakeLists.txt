INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR} ${PYTHON_INCLUDE_PATH} ${PYTHON_NUMPY_INCLUDE_DIR})

IF (CMAKE_GENERATOR MATCHES "Visual Studio")
  ADD_DEFINITIONS(-DVIGRA_DLL)
ENDIF (CMAKE_GENERATOR MATCHES "Visual Studio")

IF(HDF5_FOUND)
    ADD_DEFINITIONS(-DHasHDF5 -D_HDF5USEDLL_ -DHDF5CPP_USEDLL)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
ENDIF(HDF5_FOUND)

ADD_LIBRARY(vigranumpycmodule SHARED
    vigranumpycmodule.cxx
    vigranumpyanalysis.cxx
    vigranumpyimpex.cxx
    vigranumpykernel.cxx
    vigranumpyconverters.cxx
    vigranumpymultianalysis.cxx
    vigranumpyconvolution.cxx
    vigranumpymulticonvolution.cxx
    vigranumpyfilters.cxx
    vigranumpyrf.cxx)

IF(PYTHON_PLATFORM MATCHES "^win32$")
    SET_TARGET_PROPERTIES(vigranumpycmodule PROPERTIES PREFIX "" SUFFIX  ".pyd")
ENDIF()

TARGET_LINK_LIBRARIES(vigranumpycmodule ${Boost_PYTHON_LIBRARIES} ${PYTHON_LIBRARY} vigraimpex)
ADD_DEPENDENCIES(vigranumpy vigranumpycmodule)

IF(HDF5_FOUND)
    TARGET_LINK_LIBRARIES(vigranumpycmodule ${HDF5_LIBRARIES})
ENDIF(HDF5_FOUND)

INSTALL(TARGETS vigranumpycmodule
        RUNTIME DESTINATION ${VIGRA_NUMPY_INSTALL_DIR})

IF(FFTW3_FOUND)
    INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})
    ADD_LIBRARY(vigranumpyfourier SHARED
                vigranumpyfourier.cxx)
    TARGET_LINK_LIBRARIES(vigranumpyfourier ${FFTW3_LIBRARIES} ${Boost_PYTHON_LIBRARIES} ${PYTHON_LIBRARY})
    ADD_DEPENDENCIES(vigranumpy vigranumpyfourier)
    
    IF(PYTHON_PLATFORM MATCHES "^win32$")
        SET_TARGET_PROPERTIES(vigranumpyfourier PROPERTIES PREFIX "" SUFFIX  ".pyd")
    ENDIF()
    INSTALL(TARGETS vigranumpyfourier
            RUNTIME DESTINATION ${VIGRA_NUMPY_INSTALL_DIR})
ENDIF(FFTW3_FOUND)
