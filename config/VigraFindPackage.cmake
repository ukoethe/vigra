# - 
#
MACRO(VIGRA_FIND_PACKAGE package)
    foreach(name ${ARGN})
        SET(${package}_NAMES ${package}_NAMES ${name})
    endforeach(name)

    FIND_PACKAGE(${package} QUIET)
    
    foreach(path ${DEPENDENCY_SEARCH_PREFIX})
        if(NOT ${${package}_FOUND} MATCHES TRUE)
            SET(CMAKE_INCLUDE_PATH ${path}/include)
            SET(CMAKE_LIBRARY_PATH ${path}/lib)
            FIND_PACKAGE(${package} QUIET)
            SET(CMAKE_INCLUDE_PATH "")
            SET(CMAKE_LIBRARY_PATH "")
        endif()
    endforeach(path)
    
    if(NOT ${${package}_FOUND} MATCHES TRUE)
        MESSAGE(STATUS "package ${package} could not be found (or incompatible version)")
    endif()
    
ENDMACRO(VIGRA_FIND_PACKAGE)    
