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
    ${OnnxRuntime_INCLUDE_DIR}
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
    frnn cugraphtestutil ${TORCH_LIBRARIES}
)
target_link_libraries(
  ActsPluginExaTrkX
  PUBLIC
  ${OnnxRuntime_LIBRARY}
)

install(
  TARGETS ActsPluginExaTrkX
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(
    DIRECTORY include/Acts
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})