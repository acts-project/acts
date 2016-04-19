cmake_minimum_required (VERSION 2.6)
project (ACTS)

# set the current version of ACTS
set (ACTS_VERSION 0.0.1)

# default installation direction
set (PREFIX "/opt/ACTS/${ACTS_VERSION}" CACHE STRING "installation directory")

# check for Eigen
if (NOT DEFINED EIGEN_INCLUDE_DIR)
  message (FATAL_ERROR "'EIGEN_INCLUDE_DIR' variable not set
  include directory of Eigen library is not set
  set it using 'cmake <path/to/ACTS> -DEIGEN_INCLUDE_DIR=path/to/Eigen")
endif()

# set installation prefix
set (CMAKE_INSTALL_PREFIX ${PREFIX})

# setup include directories
include_directories (${EIGEN_INCLUDE_DIR})
include_directories (${PROJECT_SOURCE_DIR}/include)

# specify compile flags
add_definitions (-std=c++14)
add_definitions (-DACTS_PARAMETER_DEFINITIONS_PLUGIN=\"ACTS/Utilities/detail/DefaultParameterDefinitions.h\")

# get list of all source files
file (GLOB src_files "src/*/*.cxx")
add_library(ACTSCore SHARED ${src_files})

add_subdirectory (bin)
add_subdirectory (doc)

install (TARGETS ACTSCore EXPORT ACTSTargets
	 LIBRARY DESTINATION lib
	 ARCHIVE DESTINATION lib
	 RUNTIME DESTINATION bin
	 INCLUDES DESTINATION include)

install (DIRECTORY include/ACTS DESTINATION include)

include(GenerateExportHeader)
generate_export_header (ACTSCore)
set_property (TARGET ACTSCore PROPERTY VERSION ${ACTS_VERSION})

set(ConfigPackageLocation lib/cmake/ACTS)
include(CMakePackageConfigHelpers)
write_basic_package_version_file ("${ConfigPackageLocation}/ACTSConfigVersion.cmake" VERSION ${ACTS_VERSION} COMPATIBILITY AnyNewerVersion)
configure_file(cmake/ACTSConfig.cmake "${ConfigPackageLocation}/ACTSConfig.cmake" COPYONLY)
install(EXPORT ACTSTargets FILE ACTSTargets.cmake DESTINATION ${ConfigPackageLocation})

