## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
##
## # The following are required to submit to the CDash dashboard:
##   ENABLE_TESTING()
##   INCLUDE(CTest)

set(CTEST_PROJECT_NAME "Acts")
set(CTEST_NIGHTLY_START_TIME "00:00:00 EST")

# CDash configuration
set (CTEST_DROP_METHOD "http")
set (CTEST_DROP_SITE "my.cdash.org")
set (CTEST_DROP_LOCATION "/submit.php?project=Acts")
set (CTEST_DROP_SITE_CDASH TRUE)

# sub projects
set (CTEST_PROJECT_SUBPROJECTS ACore Examples DD4hepPlugin MaterialPlugin TGeoPlugin)

# do not count unit tests for test coverage
set(CTEST_CUSTOM_COVERAGE_EXCLUDE "Tests")

# site name
execute_process (COMMAND hostname OUTPUT_VARIABLE HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process (COMMAND whoami OUTPUT_VARIABLE USER OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process (COMMAND printf %s@%s ${USER} ${HOSTNAME} OUTPUT_VARIABLE CTEST_SITE)

