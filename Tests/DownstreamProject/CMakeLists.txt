cmake_minimum_required(VERSION 3.11)

project(ActsDownstreamProject)
enable_testing()

# find all optional components that are build
find_package(
    Acts
    CONFIG
    REQUIRED
    COMPONENTS Core Fatras PluginJson PluginRoot PluginCovfie
)

# place artifacts in GNU-like paths, e.g. binaries in `<build>/bin`
include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}"
)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}"
)

# link with all optional components even when they are not really used
# to check e.g. for possible linker issues
add_executable(ShowActsVersion ShowActsVersion.cpp)
target_link_libraries(
    ShowActsVersion
    PRIVATE
        Acts::Core
        Acts::Fatras
        Acts::PluginJson
        Acts::PluginRoot
        Acts::PluginCovfie
)

option(DD4HEP "Build with DD4hep" ON)
if(DD4HEP)
    message(STATUS "Adding DD4hep plugin")
    find_package(Acts CONFIG REQUIRED COMPONENTS PluginDD4hep)
    target_link_libraries(ShowActsVersion PRIVATE Acts::PluginDD4hep)
endif()

option(EDM4HEP "Build with EDM4hep" ON)
if(EDM4HEP)
    message(STATUS "Adding EDM4hep plugin")
    find_package(Acts CONFIG REQUIRED COMPONENTS PluginEDM4hep)
    target_link_libraries(ShowActsVersion PRIVATE Acts::PluginEDM4hep)
endif()

option(GEOMODEL "Build with GeoModel" ON)
if(GEOMODEL)
    message(STATUS "Adding GeoModel plugin")
    find_package(Acts CONFIG REQUIRED COMPONENTS PluginGeoModel)
    target_link_libraries(ShowActsVersion PRIVATE Acts::PluginGeoModel)
endif()

find_package(Boost REQUIRED COMPONENTS program_options)

include(ActsGeometryModuleHelpers)

acts_add_geometry_module(SimpleGeometryModule GeometryModule.cpp)

add_executable(LoadGeometryModule LoadGeometryModule.cpp)
target_link_libraries(
    LoadGeometryModule
    PRIVATE Acts::Core Boost::program_options
)

add_test(NAME DownstreamShowActsVersion COMMAND ShowActsVersion)
add_test(
    NAME DownstreamLoadGeometryModule
    COMMAND LoadGeometryModule $<TARGET_FILE:SimpleGeometryModule>
)

if(DD4HEP)
    include(ActsDD4hepGeometryModuleHelpers)

    acts_add_dd4hep_geometry_module(
        DD4hepGeometryModule
        DD4hepGeometryModule.cpp
    )

    target_link_libraries(LoadGeometryModule PRIVATE Acts::PluginDD4hep)
    target_compile_definitions(LoadGeometryModule PRIVATE ACTS_HAVE_DD4HEP)
endif()
