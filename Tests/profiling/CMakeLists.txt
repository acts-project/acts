# define executables
add_executable (profile_gsl_stepper gsl_stepper.cpp)
target_link_libraries (profile_gsl_stepper PRIVATE ACTSCore)
target_link_libraries (profile_gsl_stepper PRIVATE gsl gslcblas)

add_executable (profile_atlas_stepper atlas_stepper.cpp)
target_link_libraries (profile_atlas_stepper PRIVATE ACTSCore)

add_executable (profile_atlas_propagator atlas_propagator.cpp)
target_link_libraries (profile_atlas_propagator PRIVATE ACTSCore)

install (TARGETS profile_gsl_stepper profile_atlas_stepper profile_atlas_propagator
	     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/profiling)
