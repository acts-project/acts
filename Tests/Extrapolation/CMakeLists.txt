add_executable (MPL_checks MPL_checks)
target_link_libraries (MPL_checks PRIVATE ACTSCore)
add_test(NAME mpl_unit_test COMMAND MPL_check)

add_executable (propagator_validation propagator_validation)
target_link_libraries (propagator_validation PRIVATE ACTSCore)
target_link_libraries (propagator_validation PRIVATE gsl gslcblas)
add_test(NAME propagator_validation COMMAND propagator_validation)