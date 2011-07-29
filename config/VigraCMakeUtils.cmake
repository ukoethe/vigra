# - find the path of a module's compilation output (usually a DLL)
#
#  DEPENDENCY_PATH(TARGET_VAR dependency)
#
# e.g. DEPENDENCY_PATH(VIGRAIMPEX_PATH vigraimpex)
#
MACRO(DEPENDENCY_PATH variable target)
    GET_TARGET_PROPERTY(${variable} ${target} LOCATION)
    STRING(REGEX REPLACE "(/\\$\\([^\\)]*\\)/[^/]*|/[^/]*)$" "" ${variable} ${${variable}}) # get path prefix
    FILE(TO_NATIVE_PATH ${${variable}} ${variable})
ENDMACRO(DEPENDENCY_PATH)

