
#       message(status "Current Source: ${CMAKE_CURRENT_SOURCE_DIR}") 
#       message(status "Current Bin: ${CMAKE_CURRENT_BINARY_DIR}") 
#       configure_file(${CMAKE_CURRENT_SOURCE_DIR}/oldClassifier.log 
#                      ${CMAKE_CURRENT_BINARY_DIR}/oldClassifier.log
#                      COPYONLY)
VIGRA_COPY_TEST_DATA(oldClassifier.log oldsetTest.log)
