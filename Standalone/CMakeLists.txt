cmake_minimum_required (VERSION 2.6)
project (ACTS)

# set the current version of ACTS
set (VERSION "0.0.1")

# default installation direction
set (PREFIX "/opt/ACTS/${VERSION}" CACHE STRING "installation directory")

# check for Eigen
if (NOT DEFINED EIGEN_INC)
  message (FATAL_ERROR "'EIGEN_INC' variable not set
  include directory of Eigen library is not set
  set it using 'cmake <path/to/ACTS> -DEIGEN_INC=path/to/Eigen")
endif()

# set installation prefix
set (CMAKE_INSTALL_PREFIX ${PREFIX})

# setup include directories
include_directories (${EIGEN_INC})
include_directories (${PROJECT_SOURCE_DIR}/include)

# specify compile flags
add_definitions (-std=c++14)
add_definitions (-DACTS_PARAMETER_DEFINITIONS_PLUGIN=\"Core/DefaultParameterDefinitions.h\")

# get list of all source files
file (GLOB src_files "src/*/*.cxx")
add_library(ACTSCore SHARED ${src_files})

file(GLOB DEPLOY_FILES_AND_DIRS "${PROJECT_SOURCE_DIR}/include/*")
foreach(ITEM ${DEPLOY_FILES_AND_DIRS})
  IF( IS_DIRECTORY "${ITEM}" )
    LIST( APPEND DIRS_TO_DEPLOY "${ITEM}" )
  ELSE()
    LIST( APPEND FILES_TO_DEPLOY "${ITEM}" )
  ENDIF()
endforeach()

install (FILES ${FILES_TO_DEPLOY} DESTINATION ACTS)
install (DIRECTORY ${DIRS_TO_DEPLOY} DESTINATION ACTS)
install (TARGETS ACTSCore DESTINATION lib)
