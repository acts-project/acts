add_subdirectory(libFRNN)
add_subdirectory(cuGraph)

add_library(
  ActsPluginExaTrkX SHARED
  src/ExaTrkXTrackFinding.cpp)

# set_target_properties(ActsPluginExaTrkX PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_include_directories(
  ActsPluginExaTrkX
  PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    ${OnnxRuntime_INCLUDE_DIR}
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(
  ActsPluginExaTrkX
  PUBLIC
    frnn cugraphtestutil ${TORCH_LIBRARIES} ${OnnxRuntime_LIBRARY}
)

install(
  TARGETS ActsPluginExaTrkX
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(
    DIRECTORY include/Acts
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})