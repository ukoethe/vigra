VIGRA_ADD_TEST(test_objectfeatures test.cxx)
IF(WITH_LEMON)
    VIGRA_ADD_TEST(test_objectfeatures_lemon test_lemon.cxx LIBRARIES ${LEMON_LIBRARY})
    INCLUDE_DIRECTORIES(${LEMON_INCLUDE_DIR})
    SET_TARGET_PROPERTIES(test_objectfeatures_lemon PROPERTIES COMPILE_FLAGS "-DWITH_LEMON")
ENDIF(WITH_LEMON)
VIGRA_ADD_TEST(test_stand_alone_acc_chain stand_alone_acc_chain.cxx)
VIGRA_COPY_TEST_DATA(of.gif)

