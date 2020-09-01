# This module needs at least CMake 3.13.
cmake_minimum_required( VERSION 3.13 )

find_package(SYCL REQUIRED)

add_library(
  ActsPluginSycl SHARED
  # header files
  include/Acts/Plugins/Sycl/Seeding/detail/Types.h
  include/Acts/Plugins/Sycl/Seeding/CreateSeedsForGroupSycl.h
  include/Acts/Plugins/Sycl/Seeding/DeviceExperimentCuts.hpp
  include/Acts/Plugins/Sycl/Seeding/Seedfinder.hpp
  include/Acts/Plugins/Sycl/Seeding/Seedfinder.ipp
  include/Acts/Plugins/Sycl/Utilities/CreateQueue.h
  include/Acts/Plugins/Sycl/Utilities/DeviceSelector.h
  include/Acts/Plugins/Sycl/Utilities/ListPlatforms.h
  # implementation files
  src/Seeding/CreateSeedsForGroupSycl.cpp
  src/Utilities/CreateQueue.cpp
  src/Utilities/ListPlatforms.cpp
  src/Utilities/DeviceSelector.cpp
)

target_include_directories(
  ActsPluginSycl
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  
target_link_libraries(
  ActsPluginSycl
  PUBLIC ActsCore)

acts_setup_sycl_target(ActsPluginSycl DEPENDENCY PRIVATE)

install(
  TARGETS ActsPluginSycl
  EXPORT ActsPluginSyclTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  