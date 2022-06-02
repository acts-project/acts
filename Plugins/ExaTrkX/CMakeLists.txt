#add_subdirectory(libFRNN)

add_library(
  ActsPluginExaTrkX SHARED
  src/ExaTrkXUtils.cpp
  src/ExaTrkXTrackFindingBase.cpp
  src/ExaTrkXTrackFinding.cpp
  #src/ExaTrkXTriton.cpp
  #src/ExaTrkXTrackFindingTritonTorch.cpp
  #src/ExaTrkXTrackFindingTritonPython.cpp
  #src/ExaTrkXTrackFindingTritonOne.cpp
  #src/ExaTrkXTrackFindingTriton.cpp
)

set_target_properties(ActsPluginExaTrkX
PROPERTIES  BUILD_RPATH    "\$ORIGIN"
            INSTALL_RPATH  "\$ORIGIN"
            CXX_STANDARD                        17
            CXX_STANDARD_REQUIRED               ON
            CUDA_STANDARD                       17
            CUDA_STANDARD_REQUIRED              ON
            INTERFACE_POSITION_INDEPENDENT_CODE ON
            CUDA_SEPARABLE_COMPILATION          ON
)
target_compile_options(ActsPluginExaTrkX
    PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:${EXATRKX_CXX_FLAGS}>"
            "$<$<COMPILE_LANGUAGE:CUDA>:${EXATRKX_CUDA_FLAGS}>")
target_compile_definitions(ActsPluginExaTrkX PUBLIC CUDA_API_PER_THREAD_DEFAULT_STREAM)
target_compile_definitions(ActsPluginExaTrkX PUBLIC TRITON_ENABLE_GPU)
target_include_directories(
  ActsPluginExaTrkX
  PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_include_directories(
  ActsPluginExaTrkX
  PUBLIC
     #/workspace/install/include
)

target_link_libraries(
  ActsPluginExaTrkX
  PUBLIC
    ${TORCH_LIBRARIES}
    Boost::boost
)

target_link_libraries(
  ActsPluginExaTrkX
  PRIVATE
    frnn  
    #/workspace/install/lib/libgrpcclient.so
    #/workspace/build/third-party/protobuf/lib/libprotobuf.a
)

install(
  TARGETS ActsPluginExaTrkX
  EXPORT ActsPluginExaTrkXTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
