set(CTEST_PROJECT_NAME "Vigra")
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")


set(CTEST_DROP_METHOD "http")

set(CTEST_DROP_SITE "rubens")
set(CTEST_DROP_LOCATION "/cdash/submit.php?project=Vigra")
set(CTEST_DROP_SITE_CDASH TRUE)


SET(MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes")
