SET(TARGETS 
    convert
    subimage
    invert
    invert_explicitly
    resize
    smooth
    palette
    profile
    pyramid
    edge
    boundarytensor
    watershed
    weightedWatersheds
    voronoi
    nnlsq
    total_variation)

ADD_CUSTOM_TARGET(examples)

FOREACH(TARGET ${TARGETS})
    ADD_EXECUTABLE(example_${TARGET} EXCLUDE_FROM_ALL ${TARGET}.cxx)
    ADD_DEPENDENCIES(examples example_${TARGET})
    TARGET_LINK_LIBRARIES(example_${TARGET} vigraimpex)
ENDFOREACH(TARGET)

FIND_PACKAGE(Qt4)
IF(QT4_FOUND)
    INCLUDE_DIRECTORIES(${QT_INCLUDE_DIR} ${QT_QTCORE_INCLUDE_DIR} ${QT_QTGUI_INCLUDE_DIR} )
    ADD_EXECUTABLE(example_multi_distance multi_distance.cxx vis_vectorial_dist.cxx)
    ADD_DEPENDENCIES(examples example_multi_distance)
    TARGET_LINK_LIBRARIES(example_multi_distance vigraimpex ${QT_QTCORE_LIBRARIES} ${QT_QTGUI_LIBRARIES})
ENDIF()

get_property(location TARGET example_multi_distance PROPERTY LOCATION)
add_custom_command(TARGET example_multi_distance
POST_BUILD        
COMMAND ${CMAKE_COMMAND} -E copy_if_different  
    ${location}
    ${CMAKE_CURRENT_SOURCE_DIR}
)


