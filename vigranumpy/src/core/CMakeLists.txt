IF(MSVC)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
ENDIF()

include(VigraConfigureThreading)
VIGRA_CONFIGURE_THREADING(REQUIRED) # We require a working threading implementation.

if(WITH_BOOST_THREAD)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    SET(VIGRANUMPY_THREAD_LIBRARIES ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${Boost_CHRONO_LIBRARY})
elseif(NOT MSVC)
    IF (CMAKE_COMPILER_IS_GNUCXX AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.7.0")
        SET(CMAKE_CXX_FLAGS "-pthread -std=c++0x ${CMAKE_CXX_FLAGS}")
    elseif(CMAKE_COMPILER_IS_GNUCXX)
        SET(CMAKE_CXX_FLAGS "-pthread -std=c++11 ${CMAKE_CXX_FLAGS}")
    else()
        SET(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
    endif()
endif()

# note special treatment of target vigranumpy_core: 
# module name is automatically changed into vigranumpycore
VIGRA_ADD_NUMPY_MODULE(core 
  SOURCES 
    vigranumpycore.cxx
    converters.cxx
    axistags.cxx
    multi_array_chunked.cxx
  LIBRARIES   
    ${VIGRANUMPY_IMPEX_LIBRARIES} ${VIGRANUMPY_THREAD_LIBRARIES}
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
    non_local_mean.cxx
  LIBRARIES 
    ${VIGRANUMPY_THREAD_LIBRARIES}
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


find_package(Boost COMPONENTS system filesystem REQUIRED)

VIGRA_ADD_NUMPY_MODULE(graphs SOURCES
    graphs.cxx
    adjacencyListGraph.cxx
    gridGraphNd.cxx
    gridGraph2d.cxx
    gridGraph3d.cxx
    gridRag.cxx
    grid_graph_implicit_edge_maps.cxx
    #eccentricity.cxx
  LIBRARIES
    ${VIGRANUMPY_IMPEX_LIBRARIES}  ${VIGRANUMPY_THREAD_LIBRARIES}
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
  VIGRANUMPY)

VIGRA_ADD_NUMPY_MODULE(histogram SOURCES
    histogram.cxx
  VIGRANUMPY)

VIGRA_ADD_NUMPY_MODULE(utilities SOURCES
    utilities.cxx
  VIGRANUMPY)


VIGRA_ADD_NUMPY_MODULE(blockwise SOURCES
    blockwise.cxx
  VIGRANUMPY)

