IF(MSVC)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
ENDIF()

# note special treatment of target vigranumpy_core: 
# module name is automatically changed into vigranumpycore
VIGRA_ADD_NUMPY_MODULE(core 
  SOURCES 
    vigranumpycore.cxx
    converters.cxx
    axistags.cxx
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
    filters.cxx
    tensors.cxx
    morphology.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(analysis SOURCES
    segmentation.cxx
    edgedetection.cxx
    interestpoints.cxx
    accumulator.cxx
    accumulator-region-singleband.cxx
    accumulator-region-multiband.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(learning SOURCES
    random_forest_old.cxx
    random_forest.cxx
    learning.cxx
  LIBRARIES   
    ${VIGRANUMPY_IMPEX_LIBRARIES}
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(colors SOURCES
    colors.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(noise SOURCES
    noise.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(geometry SOURCES
    geometry.cxx
  VIGRANUMPY)
   
VIGRA_ADD_NUMPY_MODULE(optimization SOURCES
    optimization.cxx
  VIGRANUMPY)
