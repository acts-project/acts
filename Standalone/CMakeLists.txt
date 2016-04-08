cmake_minimum_required (VERSION 2.6)
project (ACTS)

LIST(APPEND CMAKE_MODULE_PATH "/var/data/a-tracking-sw/Standalone/")
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})

find_package (Eigen3 REQUIRED)

add_definitions ( -std=c++14 )
include_directories ( ${PROJECT_SOURCE_DIR}/include ${EIGEN3_INCLUDE_DIR} )

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

INSTALL( FILES ${FILES_TO_DEPLOY} DESTINATION ACTS )
INSTALL( DIRECTORY ${DIRS_TO_DEPLOY} DESTINATION ACTS )
install (TARGETS ACTSCore DESTINATION lib)
