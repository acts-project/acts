# define executables
add_executable (profile_gsl_stepper gsl_stepper.cpp)
target_link_libraries (profile_gsl_stepper PRIVATE ACTSCore)
target_link_libraries (profile_gsl_stepper PRIVATE gsl gslcblas)

install (TARGETS profile_gsl_stepper
	     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/profiling)
