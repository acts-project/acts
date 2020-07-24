
# CUDA code handling required at least CMake 3.10. Note that at the moment the
# plugin does a bad job with hiding its CUDA implementation from the outside.
# For this reason CMake needs to be told to put the CUDA Runtime API as a public
# dependency on the shared library. This will be removed after some cleanup.
cmake_minimum_required( VERSION 3.10 )

add_library(
  ActsCudaPlugin SHARED
  src/Seeding/Kernels.cu
  src/Utilities/Arrays.cu
  src/Utilities/ErrorCheck.cu
  )

target_include_directories(
  ActsCudaPlugin
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# The unit test needs to explicitly use the CUDA Runtime API.
find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})

target_link_libraries(
  ActsCudaPlugin
  PUBLIC ActsCore ${CUDART_LIBRARY})

install(
  TARGETS ActsCudaPlugin
  EXPORT ActsCudaPluginTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
