# we need boost
find_package (Boost ${ACTS_BOOST_VERSION} REQUIRED COMPONENTS program_options)

# define executables
add_executable (AtlasStepperProfiles AtlasStepperProfiles.cpp)
target_include_directories (AtlasStepperProfiles PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (AtlasStepperProfiles PRIVATE ACTSCore)
target_link_libraries (AtlasStepperProfiles PRIVATE ${Boost_LIBRARIES})

add_executable (EigenStepperProfiles EigenStepperProfiles.cpp)
target_include_directories (EigenStepperProfiles PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (EigenStepperProfiles PRIVATE ACTSCore)
target_link_libraries (EigenStepperProfiles PRIVATE ${Boost_LIBRARIES})

add_executable (AtlasPropagatorProfiles AtlasPropagatorProfiles.cpp)
target_include_directories (AtlasPropagatorProfiles PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries (AtlasPropagatorProfiles PRIVATE ACTSCore)
target_link_libraries (AtlasPropagatorProfiles PRIVATE ${Boost_LIBRARIES})

install (TARGETS AtlasStepperProfiles AtlasPropagatorProfiles
	     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/profiling)
