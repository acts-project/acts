cmake_minimum_required(VERSION 3.11)

project(ActsDownstreamProject)

# find all optional components that are build with ACTS_BUILD_EVERYTHING=on
find_package(
  Acts CONFIG REQUIRED
  COMPONENTS
    Core
    Fatras
    PluginAutodiff
    PluginDD4hep
    PluginJson
    PluginLegacy
    PluginTGeo)

# place artifacts in GNU-like paths, e.g. binaries in `<build_dir>/bin`
include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# link with all optional components even when they are not really used
# to check e.g. for possible linker issues
add_executable(ShowActsVersion ShowActsVersion.cpp)
target_link_libraries(
  ShowActsVersion
  PRIVATE
    ActsCore
    ActsFatras
    ActsPluginDD4hep
    ActsPluginJson
    ActsPluginLegacy
    ActsPluginTGeo)
