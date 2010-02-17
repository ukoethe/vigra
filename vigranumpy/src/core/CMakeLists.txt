INCLUDE_DIRECTORIES(${VIGRANUMPY_INCLUDE_DIRS})

IF(HDF5_FOUND)
    ADD_DEFINITIONS(${HDF5_CPPFLAGS})
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
    SET(VIGRANUMPY_CORE_LIBRARIES ${VIGRANUMPY_LIBRARIES} vigraimpex ${HDF5_LIBRARIES})
ELSE()
    SET(VIGRANUMPY_CORE_LIBRARIES ${VIGRANUMPY_LIBRARIES} vigraimpex)
ENDIF(HDF5_FOUND)

VIGRA_ADD_NUMPY_MODULE(core 
  SOURCES 
    vigranumpycore.cxx
    vigranumpyimagetransformation.cxx
    vigranumpyconverters.cxx
  LIBRARIES
    ${VIGRANUMPY_CORE_LIBRARIES})

VIGRA_ADD_NUMPY_MODULE(impex 
  SOURCES
    vigranumpyimpex.cxx
  LIBRARIES   
    vigraimpex ${VIGRANUMPY_LIBRARIES} vigranumpy_core)   
     
VIGRA_ADD_NUMPY_MODULE(convolution 
  SOURCES
    vigranumpykernel.cxx
    vigranumpyconvolution.cxx
    vigranumpymulticonvolution.cxx)
   
VIGRA_ADD_NUMPY_MODULE(tensor SOURCES
    vigranumpytensor.cxx)
   
VIGRA_ADD_NUMPY_MODULE(filters SOURCES
    vigranumpyfilters.cxx)
   
VIGRA_ADD_NUMPY_MODULE(noise SOURCES
    vigranumpynoise.cxx)
   
VIGRA_ADD_NUMPY_MODULE(segmentation SOURCES
    vigranumpysegmentation.cxx)
   
VIGRA_ADD_NUMPY_MODULE(morphology SOURCES
    vigranumpymorphology.cxx)
   
VIGRA_ADD_NUMPY_MODULE(edgedetection SOURCES
    vigranumpyedgedetection.cxx)
   
VIGRA_ADD_NUMPY_MODULE(classification SOURCES
    vigranumpyrf.cxx
    vigranumpyrf_new.cxx)

