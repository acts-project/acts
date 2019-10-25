configure_file(
  ActsVersion.hpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/Acts/ActsVersion.hpp)

file(GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp" "include/*.ipp")

add_library(ActsCore SHARED ${src_files})
target_compile_features(ActsCore PUBLIC cxx_std_17)
target_compile_definitions(
  ActsCore
  PUBLIC -DACTS_PARAMETER_DEFINITIONS_PLUGIN="${ACTS_PARAMETER_DEFINITIONS_PLUGIN}")
target_include_directories(
  ActsCore
  SYSTEM PUBLIC ${Boost_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS}
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         # for the generated version header
         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(ActsCore PUBLIC Boost::boost)

install(
  TARGETS ActsCore
  EXPORT ActsCoreTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/Acts/ActsVersion.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Acts)
