file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp" "include/*.ipp")

add_library(ActsCore SHARED ${src_files})

target_include_directories (ActsCore SYSTEM PUBLIC ${Boost_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS})
target_include_directories (ActsCore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_compile_definitions (ActsCore PUBLIC -DACTS_PARAMETER_DEFINITIONS_PLUGIN="${ACTS_PARAMETER_DEFINITIONS_PLUGIN}")

install (TARGETS ActsCore
         EXPORT ActsCoreTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
         COMPONENT ActsCore)

install (DIRECTORY include/Acts
         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
         COMPONENT Core)

# generate version header and make it available to the build
configure_file(
  ActsVersion.hpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/Acts/ActsVersion.hpp)

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/Acts/ActsVersion.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Acts)

target_include_directories(ActsCore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

acts_add_targets_to_cdash_project(PROJECT ACore TARGETS ActsCore)
