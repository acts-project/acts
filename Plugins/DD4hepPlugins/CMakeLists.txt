# check for Eigen
if (NOT DEFINED EIGEN_INCLUDE_DIR)
message (FATAL_ERROR "'EIGEN_INCLUDE_DIR' variable not set
include directory of Eigen library is not set
set it using 'cmake <path/to/ACTS> -DEIGEN_INCLUDE_DIR=path/to/Eigen")
endif()
# check for DD4hep
if (NOT DEFINED DD4HEP_INCLUDE_DIR)
message (FATAL_ERROR "'DD4HEP_INCLUDE_DIR' variable not set
set it using 'cmake <path/to/ACTS> -DDD4HEP_INCLUDE_DIR'")
endif()

#find_package(DD4hep)

# setup include directories
include_directories (${EIGEN_INCLUDE_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories (${DD4HEP_INCLUDE_DIR}/DDCore/include)
include_directories (${DD4HEP_INCLUDE_DIR}/DDSegmentation/include)
include_directories (${PROJECT_SOURCE_DIR}/Core/include)
include_directories (${PROJECT_SOURCE_DIR}/Plugins/TGeoPlugins/include)


# get list of all source files
file (GLOB src_files "src/*.cxx")
add_library(ACTSDD4hepPlugins SHARED ${src_files})

#target_include_directories(ACTSDD4hepPlugins )
#target_link_libraries(ACTSDD4hepPlugins )

target_include_directories(ACTSDD4hepPlugins PUBLIC Plugins/TGeoPlugins/include)
target_link_libraries(ACTSDD4hepPlugins ACTSTGeoPlugins)

install (TARGETS ACTSDD4hepPlugins EXPORT ACTSTargets
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION bin
INCLUDES DESTINATION include)

install (DIRECTORY include/ACTS DESTINATION include)

