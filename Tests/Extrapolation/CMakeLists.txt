add_executable (gsl_test gsl_test.cpp)
target_link_libraries (gsl_test PRIVATE gsl gslcblas)
add_test(NAME gsl_unit_test COMMAND gsl_test)
