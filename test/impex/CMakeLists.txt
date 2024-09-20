IF(JPEG_FOUND)
  ADD_DEFINITIONS(-DHasJPEG)
ENDIF(JPEG_FOUND)

IF(PNG_FOUND)
  ADD_DEFINITIONS(-DHasPNG)
ENDIF(PNG_FOUND)

IF(TIFF_FOUND)
  ADD_DEFINITIONS(-DHasTIFF)
  INCLUDE_DIRECTORIES(${SUPPRESS_WARNINGS} ${TIFF_INCLUDE_DIR})
ENDIF(TIFF_FOUND)

IF(OPENEXR_FOUND)
  ADD_DEFINITIONS(-DHasEXR)
ENDIF(OPENEXR_FOUND)

VIGRA_ADD_TEST(test_impex test.cxx LIBRARIES vigraimpex)

VIGRA_COPY_TEST_DATA(astronaut.xv astronaut_gifref.xv astronautfloat.xv astronautfloatrgb.xv astronautrgb.xv no-image.txt astronaut_0.tif astronaut_1.tif astronaut_2.tif astronaut_masked_color.tif  astronaut_masked_gray.tif bilevel.tiff)

