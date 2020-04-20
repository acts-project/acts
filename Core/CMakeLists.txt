configure_file(
  ActsVersion.hpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/Acts/ActsVersion.hpp)

# source files will be added later
add_library(
  ActsCore SHARED "")
target_compile_features(
  ActsCore
  PUBLIC cxx_std_17)
target_include_directories(
  ActsCore
  SYSTEM PUBLIC ${Boost_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS}
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    # for the generated version header
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(
  ActsCore
  PUBLIC Boost::boost)

if(ACTS_PARAMETER_DEFINITIONS_HEADER)
  target_compile_definitions(
    ActsCore
    PUBLIC -DACTS_PARAMETER_DEFINITIONS_HEADER="${ACTS_PARAMETER_DEFINITIONS_HEADER}")
endif()

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

# target source files are added separately
add_subdirectory(src/EventData)
add_subdirectory(src/Geometry)
add_subdirectory(src/MagneticField)
add_subdirectory(src/Material)
add_subdirectory(src/Propagator)
add_subdirectory(src/Surfaces)
add_subdirectory(src/Utilities)
add_subdirectory(src/Vertexing)
add_subdirectory(src/Visualization)
