add_subdirectory(libFRNN)
add_subdirectory(cuGraph)


add_library(
  ActsExamplesTrackFindingMLBased SHARED
  src/TrackFindingMLBasedAlgorithm.cpp)

# set_target_properties(ActsExamplesTrackFindingMLBased PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_include_directories(
  ActsExamplesTrackFindingMLBased
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         ${OnnxRuntime_INCLUDE_DIR})

target_link_libraries(
  ActsExamplesTrackFindingMLBased
  PUBLIC
    ActsCore
    ActsExamplesFramework
    frnn cugraphtestutil ${TORCH_LIBRARIES} ${OnnxRuntime_LIBRARY}
)

install(
  TARGETS ActsExamplesTrackFindingMLBased
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})