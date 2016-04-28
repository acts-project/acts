cmake_minimum_required (VERSION 2.6)
project (ACTS)

# set the current version of ACTS
set (ACTS_VERSION 0.0.1)

# default installation direction
set (CMAKE_INSTALL_PREFIX "/opt/ACTS/${ACTS_VERSION}" CACHE STRING "installation directory")

# default parameter definitions
set (ACTS_PARAMETER_DEFINITIONS_PLUGIN "ACTS/Utilities/detail/DefaultParameterDefinitions.h" CACHE STRING "default track parameter definition")

# default message macros
set (ACTS_MSG_MACROS_PLUGIN "ACTS/Utilities/detail/DefaultMsgMacros.h" CACHE STRING "default message macros")

# set location of cmake configuration files
set(ConfigPackageLocation lib/cmake/ACTS)

add_subdirectory (Core)
add_subdirectory (doc)
add_subdirectory (Plugins)

include(GenerateExportHeader)
generate_export_header (ACTSCore)
set_property (TARGET ACTSCore PROPERTY VERSION ${ACTS_VERSION})

include(CMakePackageConfigHelpers)
write_basic_package_version_file ("${CMAKE_INSTALL_PREFIX}/${ConfigPackageLocation}/ACTSConfigVersion.cmake" VERSION ${ACTS_VERSION} COMPATIBILITY AnyNewerVersion)
configure_file(cmake/ACTSConfig.cmake "${CMAKE_INSTALL_PREFIX}/${ConfigPackageLocation}/ACTSConfig.cmake" COPYONLY)
install(EXPORT ACTSTargets FILE ACTSTargets.cmake DESTINATION ${ConfigPackageLocation})

