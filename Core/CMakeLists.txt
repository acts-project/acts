find_package (Boost ${ACTS_BOOST_VERSION} REQUIRED)
find_package (Eigen ${ACTS_EIGEN_VERSION} REQUIRED)

# get list of all source files
file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp" "include/*.ipp")

# define library target
add_library (ACTSCore SHARED ${src_files})

if (SUPPRESS_LIBRARY_WARNINGS)
  target_include_directories (ACTSCore SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
  target_include_directories (ACTSCore SYSTEM PUBLIC ${EIGEN_INCLUDE_DIRS})
  message(STATUS "Suppressing warnings from libraries BOOST and EIGEN")
else (SUPPRESS_LIBRARY_WARNINGS)
  target_include_directories (ACTSCore PUBLIC ${Boost_INCLUDE_DIRS})
  target_include_directories (ACTSCore PUBLIC ${EIGEN_INCLUDE_DIRS})
endif (SUPPRESS_LIBRARY_WARNINGS)

target_include_directories (ACTSCore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_compile_definitions (ACTSCore PUBLIC -DACTS_PARAMETER_DEFINITIONS_PLUGIN="${ACTS_PARAMETER_DEFINITIONS_PLUGIN}")

install (TARGETS ACTSCore
         EXPORT ACTSCoreTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
         COMPONENT Core)
install (DIRECTORY include/ACTS
         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
         COMPONENT Core)

# add to CDash ACore subproject
acts_add_targets_to_cdash_project(PROJECT ACore TARGETS ACTSCore)
