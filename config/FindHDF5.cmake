# - Find HDF5, a library for reading and writing self describing array data.
#

FIND_PATH(HDF5_INCLUDE_DIR hdf5.h PATH_SUFFIXES hdf5/serial)

if(HDF5_INCLUDE_DIR)
    SET(HDF5_TRY_COMPILE_INCLUDE_DIR "-DINCLUDE_DIRECTORIES:STRING=${HDF5_INCLUDE_DIR}")

    set(HDF5_SUFFICIENT_VERSION FALSE)
    TRY_COMPILE(HDF5_SUFFICIENT_VERSION
                ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/config/checkHDF5version.c
                CMAKE_FLAGS "${HDF5_TRY_COMPILE_INCLUDE_DIR}")

    if(NOT HDF5_SUFFICIENT_VERSION)
        MESSAGE(STATUS "   HDF5: unable to compile a simple test program.\n      (include path: '${HDF5_INCLUDE_DIR}')" )
    else()
        SET(HDF5_VERSION_MAJOR 1)
        SET(HDF5_VERSION_MINOR 8)
        set(HDF5_SUFFICIENT_VERSION FALSE)
        TRY_COMPILE(HDF5_SUFFICIENT_VERSION
                    ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/config/checkHDF5version.c
                    COMPILE_DEFINITIONS "-DCHECK_VERSION=1 -DMIN_MAJOR=${HDF5_VERSION_MAJOR} -DMIN_MINOR=${HDF5_VERSION_MINOR}"
                    CMAKE_FLAGS "${HDF5_TRY_COMPILE_INCLUDE_DIR}")

        if(NOT HDF5_SUFFICIENT_VERSION)
            MESSAGE(STATUS "   HDF5: need at least version ${HDF5_VERSION_MAJOR}.${HDF5_VERSION_MINOR}" )
        else()
            MESSAGE(STATUS
                   "   Checking HDF5 version (at least ${HDF5_VERSION_MAJOR}.${HDF5_VERSION_MINOR}): ok")
        endif()
    endif()

    # Only configure HDF5 if a suitable version of the library was found
    if(HDF5_SUFFICIENT_VERSION)

        FIND_LIBRARY(HDF5_CORE_LIBRARY NAMES hdf5dll hdf5 hdf5dll_D hdf5_D PPATH_SUFFIXES hdf5/serial )
        FIND_LIBRARY(HDF5_HL_LIBRARY NAMES hdf5_hldll hdf5_hl hdf5_hldll_D hdf5_hl_D PATH_SUFFIXES hdf5/serial )

        # FIXME: as of version 1.8.9 and 1.8.10-patch1 (but NOT 1.8.10), these flags are
        #        already set correctly => remove or set conditionally according to version
        IF(WIN32 AND HDF5_CORE_LIBRARY MATCHES "dll.lib$")
            SET(HDF5_CFLAGS "-D_HDF5USEDLL_")
            SET(HDF5_CPPFLAGS "-D_HDF5USEDLL_ -DHDF5CPP_USEDLL")
        ELSE()
            SET(HDF5_CFLAGS)
            SET(HDF5_CPPFLAGS)
        ENDIF()
        SET(HDF5_CFLAGS "${HDF5_CPPFLAGS} -DH5_BUILT_AS_DYNAMIC_LIB")
        SET(HDF5_CPPFLAGS "${HDF5_CPPFLAGS} -DH5_BUILT_AS_DYNAMIC_LIB")


        set(HDF5_USES_ZLIB FALSE)
        TRY_COMPILE(HDF5_USES_ZLIB
                   ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/config/checkHDF5usesCompression.c
                   COMPILE_DEFINITIONS "-DH5_SOMETHING=H5_HAVE_FILTER_DEFLATE"
                   CMAKE_FLAGS "${HDF5_TRY_COMPILE_INCLUDE_DIR}")

        if(HDF5_USES_ZLIB)
            FIND_LIBRARY(HDF5_Z_LIBRARY NAMES zlib1 zlib z zlibd zlibd1)
            set(HDF5_ZLIB_OK ${HDF5_Z_LIBRARY})
        else()
            set(HDF5_ZLIB_OK TRUE)
            set(HDF5_Z_LIBRARY "")
        endif()

        set(HDF5_USES_SZLIB FALSE)
        TRY_COMPILE(HDF5_USES_SZLIB
                    ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/config/checkHDF5usesCompression.c
                    COMPILE_DEFINITIONS "-DH5_SOMETHING=H5_HAVE_FILTER_SZIP"
                    CMAKE_FLAGS "${HDF5_TRY_COMPILE_INCLUDE_DIR}")

        if(HDF5_USES_SZLIB)
            FIND_LIBRARY(HDF5_SZ_LIBRARY NAMES szlibdll sz szip)
            set(HDF5_SZLIB_OK ${HDF5_SZ_LIBRARY})
        else()
            set(HDF5_SZLIB_OK TRUE)
            set(HDF5_SZ_LIBRARY "")
        endif()
    endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set HDF5_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)

if(NOT HDF5_INCLUDE_DIR)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(HDF5 DEFAULT_MSG HDF5_INCLUDE_DIR)
elseif(NOT HDF5_SUFFICIENT_VERSION)
    # undo unsuccessful configuration
    set(HDF5_INCLUDE_DIR "")
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(HDF5 DEFAULT_MSG HDF5_SUFFICIENT_VERSION)
else()
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(HDF5 DEFAULT_MSG HDF5_CORE_LIBRARY
        HDF5_HL_LIBRARY HDF5_ZLIB_OK HDF5_SZLIB_OK HDF5_INCLUDE_DIR)
endif()

IF(HDF5_FOUND)
    SET(HDF5_LIBRARIES ${HDF5_CORE_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Z_LIBRARY} ${HDF5_SZ_LIBRARY})
ELSE()
    SET(HDF5_CORE_LIBRARY HDF5_CORE_LIBRARY-NOTFOUND)
    SET(HDF5_HL_LIBRARY   HDF5_HL_LIBRARY-NOTFOUND)
    SET(HDF5_Z_LIBRARY    HDF5_Z_LIBRARY-NOTFOUND)
    SET(HDF5_SZ_LIBRARY   HDF5_SZ_LIBRARY-NOTFOUND)
ENDIF(HDF5_FOUND)

