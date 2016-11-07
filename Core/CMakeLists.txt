# check for Eigen
if (NOT DEFINED EIGEN_INCLUDE_DIR)
  message (FATAL_ERROR "'EIGEN_INCLUDE_DIR' variable not set
  include directory of Eigen library is not set
  set it using 'cmake <path/to/ACTS> -DEIGEN_INCLUDE_DIR=path/to/Eigen")
endif (NOT DEFINED EIGEN_INCLUDE_DIR)

# we need boost
find_package (Boost 1.59 REQUIRED)

# get list of all source files
file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp" "include/*.ipp")

# define library target
add_library (ACore SHARED ${src_files})

# setup include directories
target_include_directories (ACore PUBLIC ${EIGEN_INCLUDE_DIR})
target_include_directories (ACore PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories (ACore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)

# set compiler flags
target_compile_definitions (ACore PUBLIC -DACTS_PARAMETER_DEFINITIONS_PLUGIN="${ACTS_PARAMETER_DEFINITIONS_PLUGIN}")

# set installation directories
install (TARGETS ACore EXPORT ACoreTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
         COMPONENT ACore)

install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
         COMPONENT ACore)

# add to list of supported components
list (APPEND _supported_components ACore)
set (_supported_components "${_supported_components}" PARENT_SCOPE)

add_custom_target(ACoreProject)
add_dependencies(ACoreProject ACore)