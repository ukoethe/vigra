IF(HDF5_FOUND)
    SET(VIGRANUMPY_IMPEX_LIBRARIES vigraimpex ${HDF5_LIBRARIES})
ELSE()
    SET(VIGRANUMPY_IMPEX_LIBRARIES vigraimpex)
ENDIF(HDF5_FOUND)

# note special treatment of target vigranumpy_core: 
# module name is automatically changed into vigranumpycore
VIGRA_ADD_NUMPY_MODULE(core 
  SOURCES 
    vigranumpycore.cxx
    converters.cxx
  VIGRANUMPY)

VIGRA_ADD_NUMPY_MODULE(impex 
  SOURCES
    impex.cxx
  LIBRARIES   
    ${VIGRANUMPY_IMPEX_LIBRARIES}
  VIGRANUMPY)   
     
VIGRA_ADD_NUMPY_MODULE(sampling 
  SOURCES
    sampling.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(filters SOURCES
    kernel.cxx
    convolution.cxx
    vigranumpymulticonvolution.cxx
    filters.cxx
    tensors.cxx
    morphology.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(analysis SOURCES
    segmentation.cxx
    vigranumpyedgedetection.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(learning SOURCES
    random_forest.cxx
    random_forest_new.cxx
  LIBRARIES   
    ${VIGRANUMPY_IMPEX_LIBRARIES}
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(noise SOURCES
    noise.cxx
  VIGRANUMPY)
