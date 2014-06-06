# - Try to find Chatty
# Once done this will define
#  CHATTY_FOUND - System has Chatty
#  CHATTY_INCLUDE_DIRS - The Chatty include directories
#  CHATTY_LIBRARIES - The libraries needed to use Chatty

find_path(CHATTY_INCLUDE_DIR chatty/logger.hxx)

find_library(CHATTY_LIBRARY NAMES chatty)

set(CHATTY_LIBRARIES ${CHATTY_LIBRARY} )
set(CHATTY_INCLUDE_DIRS ${CHATTY_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CHATTY_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Chatty  DEFAULT_MSG
                                  CHATTY_LIBRARY CHATTY_INCLUDE_DIR)

mark_as_advanced(CHATTY_INCLUDE_DIR CHATTY_LIBRARY)
