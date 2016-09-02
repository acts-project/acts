add_executable (gsl_test gsl_test.cpp)
target_link_libraries (gsl_test PRIVATE gsl gslcblas)
add_test(NAME gsl_unit_test COMMAND gsl_test)

add_executable (MPL_checks MPL_checks)
target_link_libraries (MPL_checks PRIVATE ACTSCore)
add_test(NAME mpl_unit_test COMMAND MPL_check)