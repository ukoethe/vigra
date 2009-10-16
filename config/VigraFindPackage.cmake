# - improved version of FIND_PACKAGE
# FIND_PACKAGE is dangerous when custom search prefixes are added, because these prefixes are tried for 
# includes and libraries independently. Thus, it may happen that a library is found under one prefix,
# whereas the corresponding include is found under another. Crashes of the compiled programs are
# then very likely. VIGRA_FIND_PACKAGE improves upon this by trying custom prefixes one at a time,
# so that either both libraries are found under the same prefix, or the search is marked as unsuccessful
# for this prefix.
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
            SET(${package}_INCLUDE_DIR ${package}_INCLUDE_DIR-NOTFOUND)
            SET(${package}_LIBRARIES ${package}_LIBRARIES-NOTFOUND)
            SET(${package}_LIBRARY ${package}_LIBRARY-NOTFOUND)
            FIND_PACKAGE(${package} QUIET)
            SET(CMAKE_INCLUDE_PATH "")
            SET(CMAKE_LIBRARY_PATH "")
        endif()
    endforeach(path)
    
    if(NOT ${${package}_FOUND} MATCHES TRUE)
        MESSAGE(STATUS "package ${package} could not be found (or incompatible version)")
    endif()
    
ENDMACRO(VIGRA_FIND_PACKAGE)    
