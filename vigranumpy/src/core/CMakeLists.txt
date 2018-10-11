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

IF(WITH_LEMON)
    VIGRA_ADD_NUMPY_MODULE(analysis SOURCES
        segmentation.cxx
        edgedetection.cxx
        interestpoints.cxx
        accumulator.cxx
        accumulator-region-singleband.cxx
        accumulator-region-multiband.cxx
      LIBRARIES
        ${VIGRANUMPY_THREAD_LIBRARIES}
        ${LEMON_LIBRARY}
      VIGRANUMPY)
    INCLUDE_DIRECTORIES(${LEMON_INCLUDE_DIR})
    SET_TARGET_PROPERTIES(vigranumpy_analysis PROPERTIES COMPILE_FLAGS "-DWITH_LEMON")
ELSE(WITH_LEMON)
    VIGRA_ADD_NUMPY_MODULE(analysis SOURCES
        segmentation.cxx
        edgedetection.cxx
        interestpoints.cxx
        accumulator.cxx
        accumulator-region-singleband.cxx
        accumulator-region-multiband.cxx
      LIBRARIES
        ${VIGRANUMPY_THREAD_LIBRARIES}
      VIGRANUMPY)
ENDIF(WITH_LEMON)

VIGRA_ADD_NUMPY_MODULE(learning SOURCES
    random_forest_old.cxx
    random_forest.cxx
    random_forest_3.cxx
    learning.cxx
  LIBRARIES
    ${VIGRANUMPY_IMPEX_LIBRARIES} ${VIGRANUMPY_THREAD_LIBRARIES}
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



VIGRA_ADD_NUMPY_MODULE(graphs SOURCES
    graphs.cxx
    adjacencyListGraph.cxx
    gridGraphNd.cxx
    gridGraph2d.cxx
    gridGraph3d.cxx
    grid_graph_implicit_edge_maps.cxx
    #eccentricity.cxx
  LIBRARIES
    ${VIGRANUMPY_IMPEX_LIBRARIES} ${VIGRANUMPY_THREAD_LIBRARIES}
  VIGRANUMPY)

VIGRA_ADD_NUMPY_MODULE(histogram SOURCES
    histogram.cxx
  VIGRANUMPY)

VIGRA_ADD_NUMPY_MODULE(utilities SOURCES
    utilities.cxx
  VIGRANUMPY)


VIGRA_ADD_NUMPY_MODULE(blockwise SOURCES
    blockwise.cxx
  LIBRARIES
    ${VIGRANUMPY_THREAD_LIBRARIES}
  VIGRANUMPY)

