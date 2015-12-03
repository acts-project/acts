CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#---------------------------------------------------------------
# Load macros and functions for Gaudi-based projects

find_package(GaudiProject)

#---------------------------------------------------------------

# Declare project name and version
# Syntax:
#    gaudi_project(this_project this_version
#                  USE dep_project_1 version_1 [project_2 version_2 ...]
#                  [DATA pkg1 [VERSION vers1] [pkg2 ...])
gaudi_project(a-tracking-sw HEAD
	      USE Gaudi v25r2 )

include_directories()

