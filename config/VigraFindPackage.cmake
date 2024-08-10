# - improved version of the reduced FIND_PACKAGE syntax
#
#  VIGRA_FIND_PACKAGE(package [VERSION_SPEC] [COMPONENTS component1 ... ] [NAMES name1 ...])
#
# FIND_PACKAGE is dangerous when custom search prefixes are added, because these
# prefixes are tried for includes and libraries independently. Thus, it may happen
# that a library is found under one prefix, whereas the corresponding include is
# found under another. Subsequent crashes of the compiled programs are likely.
#
# VIGRA_FIND_PACKAGE improves upon this by trying custom prefixes one at a time,
# so that either both includes and libraries are found under the same prefix,
# or the search is marked as unsuccessful for this prefix.
#
MACRO(VIGRA_FIND_PACKAGE package)
    # parse the args
    set(COMPONENTS)
    set(NAMES)
    set(VERSION_SPEC)
    set(v VERSION_SPEC)
    foreach(i ${ARGN})
        if(${i} MATCHES "^COMPONENTS$")
            set(v COMPONENTS)
        elseif(${i} MATCHES "^NAMES$")
            set(v NAMES)
        else()
            set(${v} ${${v}} ${i})
        endif()
    endforeach(i)

    SET(${package}_NAMES ${NAMES} ${${package}_NAMES})

    IF(DEFINED COMPONENTS)
        SET(COMPONENTS COMPONENTS ${COMPONENTS})
    ENDIF()

    if(VERSION_SPEC)
        set(VERSION_MESSAGE " (at least version ${VERSION_SPEC})")
    else()
        set(VERSION_MESSAGE)
    endif()

    MESSAGE(STATUS "Searching for ${package}${VERSION_MESSAGE}")

    FIND_PACKAGE(${package} ${VERSION_SPEC} ${COMPONENTS})

    # search the package in the default locations if not found
    # in the DEPENDENCY_SEARCH_PREFIX
    if(NOT ${package}_FOUND)
        MESSAGE(STATUS "   in default locations")
        FIND_PACKAGE(${package} ${VERSION_SPEC} ${COMPONENTS})
    endif()

ENDMACRO(VIGRA_FIND_PACKAGE)
