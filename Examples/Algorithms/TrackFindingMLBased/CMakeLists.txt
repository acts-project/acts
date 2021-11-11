
include(FetchContent)
FetchContent_Declare(
  rapids-cmake
  GIT_REPOSITORY https://github.com/rapidsai/rapids-cmake.git
  GIT_TAG        origin/branch-21.10
  )
FetchContent_MakeAvailable(rapids-cmake)

include(rapids-cpm)
include(rapids-cmake)
include(rapids-cuda)
include(rapids-export)
include(rapids-find)

rapids_cuda_init_architectures(CUGRAPH_TESTS)


if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA" AND
   CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 11.0)
    message(FATAL_ERROR "CUDA compiler version must be at least 11.0")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
   CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.3)
    message(FATAL_ERROR "GCC compiler must be at least 9.3")
endif()


# Remove the following archs from CMAKE_CUDA_ARCHITECTURES that
# cuhornet currently doesn't support
#
# >= 86
set(supported_archs "60" "62" "70" "72" "75" "80")
foreach( arch IN LISTS CMAKE_CUDA_ARCHITECTURES)
    string(REPLACE "-real" "" arch ${arch})
    if( arch IN_LIST supported_archs )
        list(APPEND usable_arch_values ${arch})
    endif()
endforeach()
# Make sure everything but the 'newest' arch
# is marked as `-real` so we only generate PTX for
# arch > 80
list(POP_BACK usable_arch_values latest_arch)
list(TRANSFORM usable_arch_values APPEND "-real")
list(APPEND usable_arch_values ${latest_arch})
set(CMAKE_CUDA_ARCHITECTURES ${usable_arch_values})


option(CUDA_STATIC_RUNTIME "Use CUDA static runtime" OFF)
option(CMAKE_CUDA_LINEINFO "Enable the -lineinfo option for nvcc (useful for cuda-memcheck / profiler" OFF)

###################################################################################################
# - build type ------------------------------------------------------------------------------------

# Set a default build type if none was specified
rapids_cmake_build_type(Release)

# Set the CUDA runtime CUDA targets will use
rapids_cuda_init_runtime(USE_STATIC ${CUDA_STATIC_RUNTIME})

# conda environment
rapids_cmake_support_conda_env( conda_env MODIFY_PREFIX_PATH )

###################################################################################################
# - compiler options ------------------------------------------------------------------------------

# Find the CUDAToolkit
rapids_find_package(CUDAToolkit REQUIRED
                    BUILD_EXPORT_SET cugraphTest-targets
                    INSTALL_EXPORT_SET cugraphTest-targets
                    )

# find Threads (needed by cudftestutil)
rapids_find_package(Threads REQUIRED
                    BUILD_EXPORT_SET cugraphTest-targets
                    INSTALL_EXPORT_SET cugraphTest-targets
                    )

###################################################################################################
# - find blas -------------------------------------------------------------------------------------

if(NOT DEFINED BLAS_LIBRARIES)
rapids_find_package(BLAS REQUIRED 
                    BUILD_EXPORT_SET cugraphTest-targets
                    INSTALL_EXPORT_SET cugraphTest-targets
              )
else()
  message(STATUS "Manually setting BLAS to ${BLAS_LIBRARIES}")
endif()

###################################################################################################
# - find CPM based dependencies  ------------------------------------------------------------------
rapids_cpm_init()
include(cmake/get_cuco.cmake)
include(cmake/get_nccl.cmake)

set(CUGRAPH_CXX_FLAGS "")
set(CUGRAPH_CUDA_FLAGS "")

if(CMAKE_COMPILER_IS_GNUCXX)
    list(APPEND CUGRAPH_CXX_FLAGS -Werror -Wno-error=deprecated-declarations)
endif(CMAKE_COMPILER_IS_GNUCXX)


message("-- Building for GPU_ARCHS = ${CMAKE_CUDA_ARCHITECTURES}")

list(APPEND CUGRAPH_CUDA_FLAGS --expt-extended-lambda --expt-relaxed-constexpr)
#list(APPEND CUGRAPH_CUDA_FLAGS -Werror=cross-execution-space-call -Wno-deprecated-declarations -Xptxas=--disable-warnings)
list(APPEND CUGRAPH_CUDA_FLAGS -Xcompiler=-Wall,-Wno-error=sign-compare,-Wno-error=unused-but-set-variable)
list(APPEND CUGRAPH_CUDA_FLAGS -Xfatbin=-compress-all)

# Option to enable line info in CUDA device compilation to allow introspection when profiling /
# memchecking
if (CMAKE_CUDA_LINEINFO)
    list(APPEND CUGRAPH_CUDA_FLAGS -lineinfo)
endif(CMAKE_CUDA_LINEINFO)

# Debug options
if(CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Building with debugging flags")
    list(APPEND CUGRAPH_CUDA_FLAGS -G -Xcompiler=-rdynamic)
endif(CMAKE_BUILD_TYPE MATCHES Debug)


set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_C_FLAGS "-Wno-error=format-truncation") 

add_subdirectory(libFRNN)
add_subdirectory(cuGraph)
find_package(OnnxRuntime)

add_library(
  ActsExamplesTrackFindingMLBased SHARED
  src/TrackFindingMLBasedAlgorithm.cpp)

set_target_properties(ActsExamplesTrackFindingMLBased PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_include_directories(
  ActsExamplesTrackFindingMLBased
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         ${OnnxRuntime_INCLUDE_DIR})

target_link_libraries(
  ActsExamplesTrackFindingMLBased
  PUBLIC
    ActsCore
    ActsExamplesFramework
    frnn cugraphtestutil ${TORCH_LIBRARIES} ${OnnxRuntime_LIBRARY}
)

install(
  TARGETS ActsExamplesTrackFindingMLBased
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})