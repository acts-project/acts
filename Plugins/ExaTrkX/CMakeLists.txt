find_package(NCCL REQUIRED)
find_package(cugraph REQUIRED)


# TODO this is just to get this work somehow, later check for better solutions
add_compile_definitions("SPDLOG_FMT_EXTERNAL")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")

add_subdirectory(libFRNN)
add_subdirectory(cuGraph)

add_library(
  ActsPluginExaTrkX SHARED
  src/ExaTrkXTrackFinding.cpp)
  
find_package(Torch REQUIRED)

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
    cugraphtestutil
    ${TORCH_LIBRARIES}
    NCCL::NCCL
)

target_link_libraries(
  ActsPluginExaTrkX
  PUBLIC
    OnnxRuntime  
)

install(
  TARGETS ActsPluginExaTrkX
  EXPORT ActsPluginExaTrkXTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(
    DIRECTORY include/Acts
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
