# - Find FFTW3F
# Find the single-precision FFTW3 includes and library
# This module defines
#  FFTW3F_INCLUDE_DIR, where to find fftw3.h, etc.
#  FFTW3F_LIBRARIES, the libraries needed to use single-precision FFTW3.
#  JFFTW3_FOUND, If false, do not try to use FFTW3.
# also defined, but not for general use are
#  FFTW3F_LIBRARY, where to find the single-precision FFTW3 library.

FIND_PATH(FFTW3F_INCLUDE_DIR fftw3.h)

SET(FFTW3F_NAMES ${FFTW3F_NAMES} fftw3f)
FIND_LIBRARY(FFTW3F_LIBRARY NAMES ${FFTW3F_NAMES} )

# handle the QUIETLY and REQUIRED arguments and set FFTW3F_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW3F DEFAULT_MSG FFTW3F_LIBRARY FFTW3F_INCLUDE_DIR)

IF(FFTW3F_FOUND)
  SET(FFTW3F_LIBRARIES ${FFTW3F_LIBRARY})
ENDIF(FFTW3F_FOUND)

# Deprecated declarations.
SET (NATIVE_FFTW3F_INCLUDE_PATH ${FFTW3F_INCLUDE_DIR} )
IF(FFTW3F_LIBRARY)
  GET_FILENAME_COMPONENT (NATIVE_FFTW3F_LIB_PATH ${FFTW3F_LIBRARY} PATH)
ENDIF(FFTW3F_LIBRARY)
