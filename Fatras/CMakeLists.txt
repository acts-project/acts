add_library(
  ActsFatras SHARED
  src/LandauDistribution.cpp
  src/Particle.cpp
  src/ParticleData.cpp
  src/ProcessType.cpp
  src/StandardPhysicsLists.cpp)
target_compile_features(
  ActsFatras
  PUBLIC cxx_std_17)
target_include_directories(
  ActsFatras
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(
  ActsFatras
  PUBLIC ActsCore)

install(
  TARGETS ActsFatras
  EXPORT ActsFatrasTargets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(
  DIRECTORY include/ActsFatras
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ActsFatras)
