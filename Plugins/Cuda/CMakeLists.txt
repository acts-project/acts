
# Set up the "CUDA 1" plugin.
find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})

add_library(
  ActsCudaPlugin SHARED
  src/Seeding/Kernels.cu
  )

target_include_directories(
  ActsCudaPlugin
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

target_link_libraries(
  ActsCudaPlugin
  PUBLIC ActsCore ${CUDART_LIBRARY})

# Set up the "CUDA 2" plugin.
add_library(
  ActsCuda2Plugin STATIC
  include/Acts/Plugins/Cuda/Seeding2/Details/CountDublets.hpp
  include/Acts/Plugins/Cuda/Seeding2/Details/FindDublets.hpp
  include/Acts/Plugins/Cuda/Seeding2/Details/FindTriplets.hpp
  include/Acts/Plugins/Cuda/Seeding2/Details/Types.hpp
  include/Acts/Plugins/Cuda/Seeding2/SeedFinder.hpp
  include/Acts/Plugins/Cuda/Seeding2/SeedFinder.ipp
  include/Acts/Plugins/Cuda/Seeding2/TripletFilterConfig.hpp
  include/Acts/Plugins/Cuda/Utilities/Arrays.hpp
  src/Seeding2/CountDublets.cu
  src/Seeding2/FindDublets.cu
  src/Seeding2/FindTriplets.cu
  src/Utilities/Arrays.cu
  src/Utilities/ErrorCheck.cuh
  src/Utilities/ErrorCheck.cu
  src/Utilities/MatrixMacros.hpp
)

target_include_directories(
  ActsCuda2Plugin
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(
  ActsCuda2Plugin
  PUBLIC ActsCore)

set_target_properties(ActsCuda2Plugin PROPERTIES
  CUDA_SEPARABLE_COMPILATION ON
  POSITION_INDEPENDENT_CODE ON)

# Install all CUDA plugins.
install(
  TARGETS ActsCudaPlugin ActsCuda2Plugin
  EXPORT ActsCudaPluginTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
