if(HDF5_FOUND)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  
    ADD_DEFINITIONS(-DHasHDF5 ${HDF5_CPPFLAGS})
endif()

VIGRA_ADD_TEST(test_sifImport test.cxx LIBRARIES vigraimpex)
VIGRA_COPY_TEST_DATA(testSif_4_16_30001.sif testSif_4_13_30000.sif testSif_4_6_30000.sif testSif_forBlocks_4_16_30001.sif)
