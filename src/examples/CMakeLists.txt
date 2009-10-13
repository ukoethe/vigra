SET(TARGETS 
    convert
    subimage
    invert
    invert_explicitly
    resize
    smooth
    profile
    pyramid
    edge
    watershed
    weightedWatersheds
    voronoi)

ADD_CUSTOM_TARGET(examples)

FOREACH(TARGET ${TARGETS})
    ADD_EXECUTABLE(${TARGET} ${TARGET}.cxx)
    ADD_DEPENDENCIES(examples ${TARGET})
    TARGET_LINK_LIBRARIES(${TARGET} vigraimpex)
ENDFOREACH(TARGET)
