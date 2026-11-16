cmake_minimum_required(VERSION 3.11)

project(ActsDownstreamProjectNodeps)

set(_consumer_minimum_version "${CMAKE_MINIMUM_REQUIRED_VERSION}")
find_package(Acts CONFIG REQUIRED COMPONENTS Core)
if(NOT CMAKE_MINIMUM_REQUIRED_VERSION STREQUAL _consumer_minimum_version)
    message(FATAL_ERROR "ActsConfig changed the consumer minimum CMake version")
endif()

# place artifacts in GNU-like paths, e.g. binaries in `<build>/bin`
include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}"
)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}"
)

add_executable(ShowActsVersion ShowActsVersion.cpp)
target_link_libraries(
    ShowActsVersion
    PRIVATE Acts::Core Boost::boost Eigen3::Eigen
)
