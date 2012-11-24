# - Find LEMON
# Find the native LEMON includes and library
# This module defines
#  LEMON_INCLUDE_DIR, where to find lemon/core.h, etc.
#  LEMON_LIBRARIES, the libraries needed to use LEMON.
# also defined, but not for general use are
#  LEMON_LIBRARY, where to find the LEMON library.

FIND_PATH(LEMON_INCLUDE_DIR lemon/core.h)

SET(LEMON_NAMES ${LEMON_NAMES} lemon)
FIND_LIBRARY(LEMON_LIBRARY NAMES ${LEMON_NAMES} )

# handle the QUIETLY and REQUIRED arguments and set LEMON_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LEMON DEFAULT_MSG LEMON_LIBRARY LEMON_INCLUDE_DIR)

IF(LEMON_FOUND)
  SET(LEMON_LIBRARIES ${LEMON_LIBRARY})
ENDIF(LEMON_FOUND)
