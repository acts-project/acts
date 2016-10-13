# we need boost
find_package (Boost 1.61 REQUIRED COMPONENTS program_options)

# define executables
#add_executable (profile_gsl_stepper gsl_stepper.cpp)
#target_include_directories (profile_gsl_stepper PRIVATE ${Boost_INCLUDE_DIRS})
#target_link_libraries (profile_gsl_stepper PRIVATE ACTSCore)
#target_link_libraries (profile_gsl_stepper PRIVATE gsl gslcblas)
#target_link_libraries (profile_gsl_stepper PRIVATE ${Boost_LIBRARIES})

add_executable (profile_atlas_stepper atlas_stepper.cpp)
target_include_directories (profile_atlas_stepper PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (profile_atlas_stepper PRIVATE ACTSCore)
target_link_libraries (profile_atlas_stepper PRIVATE ${Boost_LIBRARIES})

add_executable (profile_eigen_stepper eigen_stepper.cpp)
target_include_directories (profile_eigen_stepper PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (profile_eigen_stepper PRIVATE ACTSCore)
target_link_libraries (profile_eigen_stepper PRIVATE ${Boost_LIBRARIES})

add_executable (profile_atlas_propagator atlas_propagator.cpp)
target_include_directories (profile_atlas_propagator PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (profile_atlas_propagator PRIVATE ACTSCore)
target_link_libraries (profile_atlas_propagator PRIVATE ${Boost_LIBRARIES})

install (TARGETS profile_atlas_stepper profile_atlas_propagator
	     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/profiling)
