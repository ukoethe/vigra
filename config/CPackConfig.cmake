# general information
SET(CPACK_PACKAGE_VENDOR "Heidelberg Collaboratory for Image Processing")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "C++ computer vision library with emphasize on customizable algorithms and data structures"
)

# package version setup
STRING(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+" "\\1" CPACK_PACKAGE_VERSION_MAJOR "${vigra_version}")
STRING(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+" "\\1" CPACK_PACKAGE_VERSION_MINOR "${vigra_version}")
STRING(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+)" "\\1" CPACK_PACKAGE_VERSION_PATCH "${vigra_version}")

SET(CPACK_PACKAGE_INSTALL_DIRECTORY     "${PROJECT_NAME}")
SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY  "${PROJECT_NAME}")
SET(CPACK_RESOURCE_FILE_LICENSE         "${PROJECT_SOURCE_DIR}/LICENSE.txt")
SET(CPACK_RESOURCE_FILE_README          "${PROJECT_SOURCE_DIR}/README.txt")
SET(CPACK_STRIP_FILES TRUE)
SET(CPACK_SOURCE_IGNORE_FILES .hg ".#" "#.*~")
SET(CPACK_PACKAGE_CONTACT "Ullrich Koethe <ullrich.koethe@iwr.uni-heidelberg.de>")
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "")

INCLUDE (CPack)
