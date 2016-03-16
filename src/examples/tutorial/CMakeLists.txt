# SET(TARGETS
    # composite
    # dissolve
    # graph_agglomerative_clustering
    # imageExportInfo_tutorial
    # imageIO_tutorial
    # imageImportInfo_tutorial
    # invert_tutorial
    # mirror_tutorial
    # smooth_blockwise
    # smooth_convolve
    # smooth_explicitly
    # subimage_tutorial
    # transpose
    # transpose_image_tutorial)

SET(TARGETS
    smooth_blockwise
    graph_agglomerative_clustering
    )

FOREACH(TARGET ${TARGETS})
    ADD_EXECUTABLE(example_${TARGET} EXCLUDE_FROM_ALL ${TARGET})
    ADD_DEPENDENCIES(examples example_${TARGET})
    TARGET_LINK_LIBRARIES(example_${TARGET} vigraimpex)
ENDFOREACH(TARGET)
