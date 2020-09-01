#find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})   

# This module needs at least CMake 3.13.
cmake_minimum_required( VERSION 3.13 )

# We use CMake's built-in modules, used to look for SYCL capabilities.
include( CheckIncludeFileCXX )
include( CheckCXXSourceCompiles )
include( CMakeParseArguments )
set( CMAKE_REQUIRED_QUIET TRUE )

# set(SYCL_DIR ${CMAKE_CURRENT_SOURCE_DIR})
find_package(SYCL REQUIRED)

add_library(
  ActsPluginSycl SHARED
  src/Seeding/CreateSeedsForGroupSycl.cpp
  src/Utilities/CreateQueueSycl.cpp
  src/Utilities/ListSyclPlatforms.cpp
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