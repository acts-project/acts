# check for Eigen
if (NOT DEFINED EIGEN_INCLUDE_DIR)
message (FATAL_ERROR "'EIGEN_INCLUDE_DIR' variable not set
include directory of Eigen library is not set
set it using 'cmake <path/to/ACTS> -DEIGEN_INCLUDE_DIR=path/to/Eigen")
endif()

# setup include directories
include_directories (${EIGEN_INCLUDE_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories (${PROJECT_SOURCE_DIR}/Core/include)
include_directories ($ENV{ROOTSYS}/include)

# get list of all source files
file (GLOB src_files "src/*.cxx")
add_library(ACTSGenericDetectorPlugins SHARED ${src_files})

install (TARGETS ACTSGenericDetectorPlugins EXPORT ACTSTargets
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION bin
INCLUDES DESTINATION include)

install (DIRECTORY include/ACTS DESTINATION include)
