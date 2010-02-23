IF(HDF5_FOUND)
    SET(VIGRANUMPY_IMPEX_LIBRARIES vigraimpex ${HDF5_LIBRARIES})
ELSE()
    SET(VIGRANUMPY_IMPEX_LIBRARIES vigraimpex)
ENDIF(HDF5_FOUND)

VIGRA_ADD_NUMPY_MODULE(core 
  SOURCES 
    vigranumpycore.cxx
    vigranumpyimagetransformation.cxx
    vigranumpyconverters.cxx
  VIGRANUMPY)

VIGRA_ADD_NUMPY_MODULE(impex 
  SOURCES
    vigranumpyimpex.cxx
  LIBRARIES   
    ${VIGRANUMPY_IMPEX_LIBRARIES}
  VIGRANUMPY)   
     
VIGRA_ADD_NUMPY_MODULE(convolution 
  SOURCES
    vigranumpykernel.cxx
    vigranumpyconvolution.cxx
    vigranumpymulticonvolution.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(tensor SOURCES
    vigranumpytensor.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(filters SOURCES
    vigranumpyfilters.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(noise SOURCES
    vigranumpynoise.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(segmentation SOURCES
    vigranumpysegmentation.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(morphology SOURCES
    vigranumpymorphology.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(edgedetection SOURCES
    vigranumpyedgedetection.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(classification SOURCES
    vigranumpyrf.cxx
    vigranumpyrf_new.cxx
  LIBRARIES   
    ${VIGRANUMPY_IMPEX_LIBRARIES}
  VIGRANUMPY)

