cmake_minimum_required(VERSION 2.6)

INCLUDE_DIRECTORIES(../../include)

FIND_PACKAGE(JPEG)
IF(JPEG_FOUND)
  ADD_DEFINITIONS(-DHasJPEG)
ENDIF(JPEG_FOUND)

FIND_PACKAGE(PNG)
IF(PNG_FOUND)
  ADD_DEFINITIONS(-DHasPNG)
ENDIF(PNG_FOUND)

FIND_PACKAGE(TIFF)
IF(TIFF_FOUND)
  ADD_DEFINITIONS(-DHasTIFF)
ENDIF(TIFF_FOUND)

ADD_DEFINITIONS(-DVIGRA_DLL)

ADD_LIBRARY(vigraimpex SHARED
    bmp.cxx
    byteorder.cxx
    codecmanager.cxx
    gif.cxx
    hdr.cxx
    iccjpeg.c
    imageinfo.cxx
    jpeg.cxx
    png.cxx
    pnm.cxx
    rgbe.c
    sun.cxx
    tiff.cxx
    viff.cxx
    void_vector.cxx)

INSTALL(TARGETS vigraimpex
  RUNTIME DESTINATION bin 
  LIBRARY DESTINATION lib 
  ARCHIVE DESTINATION lib)
