add_library(
  ActsPluginExaTrkX SHARED
  src/ExaTrkXTrackFinding.cpp)

# set_target_properties(ActsPluginExaTrkX PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_include_directories(
  ActsPluginExaTrkX
  PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_include_directories(
  ActsPluginExaTrkX
  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

target_link_libraries(
  ActsPluginExaTrkX
  PRIVATE
    frnn
    ${TORCH_LIBRARIES}
)

target_link_libraries(
  ActsPluginExaTrkX
  PUBLIC
    Boost::boost
    OnnxRuntime
    cugraph::cugraph
)

install(
  TARGETS ActsPluginExaTrkX
  EXPORT ActsPluginExaTrkXTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(
    DIRECTORY include/Acts
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
