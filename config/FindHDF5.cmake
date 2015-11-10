# - Find HDF5, a library for reading and writing self describing array data.
#
FIND_PATH(HDF5_INCLUDE_DIR hdf5.h)

if(HDF5_INCLUDE_DIR)
    SET(HDF5_TRY_COMPILE_INCLUDE_DIR "-DINCLUDE_DIRECTORIES:STRING=${HDF5_INCLUDE_DIR}")

    FIND_LIBRARY(HDF5_CORE_LIBRARY NAMES hdf5dll hdf5  )
    FIND_LIBRARY(HDF5_HL_LIBRARY NAMES hdf5_hldll hdf5_hl  )

    # FIXME: as of version 1.8.9 and 1.8.10-patch1 (but NOT 1.8.10), these flags are
    #        already set correctly => remove or set conditionally according to version
    IF(WIN32 AND HDF5_CORE_LIBRARY MATCHES "dll.lib$")
        SET(HDF5_CFLAGS "-D_HDF5USEDLL_")
        SET(HDF5_CPPFLAGS "-D_HDF5USEDLL_ -DHDF5CPP_USEDLL")
    ELSE()
        SET(HDF5_CFLAGS)
        SET(HDF5_CPPFLAGS)
    ENDIF()

    SET(HDF5_VERSION_MAJOR 1)
    SET(HDF5_VERSION_MINOR 8)

    set(HDF5_SUFFICIENT_VERSION FALSE)
    TRY_COMPILE(HDF5_SUFFICIENT_VERSION 
                ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/config/checkHDF5version.c
                COMPILE_DEFINITIONS "-DMIN_MAJOR=${HDF5_VERSION_MAJOR} -DMIN_MINOR=${HDF5_VERSION_MINOR}"
                CMAKE_FLAGS "${HDF5_TRY_COMPILE_INCLUDE_DIR}") 
            
    if(HDF5_SUFFICIENT_VERSION)
        MESSAGE(STATUS 
               "Checking HDF5 version (at least ${HDF5_VERSION_MAJOR}.${HDF5_VERSION_MINOR}): ok")
    else()
        MESSAGE( STATUS "HDF5: need at least version ${HDF5_VERSION_MAJOR}.${HDF5_VERSION_MINOR}" )
    endif()
    
    if(HDF5_SUFFICIENT_VERSION)
        # Only test for HDF5 features if a suitable version of the library was 
        # found previously.

        set(HDF5_USES_ZLIB FALSE)
        TRY_COMPILE(HDF5_USES_ZLIB
                   ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/config/checkHDF5usesCompression.c
                   COMPILE_DEFINITIONS "-DH5_SOMETHING=H5_HAVE_FILTER_DEFLATE"
                   CMAKE_FLAGS "${HDF5_TRY_COMPILE_INCLUDE_DIR}") 

        if(HDF5_USES_ZLIB)
            FIND_LIBRARY(HDF5_Z_LIBRARY NAMES zlib1 zlib z )
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

FIND_PACKAGE_HANDLE_STANDARD_ARGS(HDF5 DEFAULT_MSG HDF5_CORE_LIBRARY 
        HDF5_HL_LIBRARY HDF5_ZLIB_OK HDF5_SZLIB_OK HDF5_INCLUDE_DIR HDF5_SUFFICIENT_VERSION)
        
IF(HDF5_FOUND)
    SET(HDF5_LIBRARIES ${HDF5_CORE_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Z_LIBRARY} ${HDF5_SZ_LIBRARY})
ELSE()
    SET(HDF5_CORE_LIBRARY HDF5_CORE_LIBRARY-NOTFOUND)
    SET(HDF5_HL_LIBRARY   HDF5_HL_LIBRARY-NOTFOUND)
    SET(HDF5_Z_LIBRARY    HDF5_Z_LIBRARY-NOTFOUND)
    SET(HDF5_SZ_LIBRARY   HDF5_SZ_LIBRARY-NOTFOUND)
ENDIF(HDF5_FOUND)
    
