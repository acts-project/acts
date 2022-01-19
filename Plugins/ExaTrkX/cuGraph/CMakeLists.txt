#=============================================================================
#
# Copyright (c) 2019-2021, NVIDIA CORPORATION.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#=============================================================================

###################################################################################################
# - common test utils -----------------------------------------------------------------------------

add_library(cugraphtestutil STATIC
            utilities/matrix_market_file_utilities.cu
            # utilities/rmat_utilities.cu
            utilities/thrust_wrapper.cu
            utilities/misc_utilities.cpp
            components/wcc_graphs.cu)

target_compile_options(cugraphtestutil
            PUBLIC "$<$<COMPILE_LANGUAGE:CXX>:${CUGRAPH_CXX_FLAGS}>"
                   "$<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:CUDA>:${CUGRAPH_CUDA_FLAGS}>>"
)

set_property(TARGET cugraphtestutil PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(cugraphtestutil
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/mmio"
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(cugraphtestutil
    PUBLIC
        #$<BUILD_INTERFACE:conda_env>
        cugraph::cugraph
        cuco::cuco
        NCCL::NCCL
)

####################

enable_language(C)
add_library(mmio STATIC mmio/mmio.c)
set_property(TARGET mmio PROPERTY LINKER_LANGUAGE C)

target_include_directories(mmio 
    PUBLIC 
        "${CMAKE_CURRENT_SOURCE_DIR}/mmio"
)

# Testing the weakly connected components
add_executable(test_connected_components test_connected_components.cpp)

set_target_properties(test_connected_components
PROPERTIES  BUILD_RPATH    "\$ORIGIN"
            INSTALL_RPATH  "\$ORIGIN"
            CXX_STANDARD                        17
            CXX_STANDARD_REQUIRED               ON
            CUDA_STANDARD                       17
            CUDA_STANDARD_REQUIRED              ON
            INTERFACE_POSITION_INDEPENDENT_CODE ON
)

target_compile_options(test_connected_components
    PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:${CUGRAPH_CXX_FLAGS}>"
            "$<$<COMPILE_LANGUAGE:CUDA>:${CUGRAPH_CUDA_FLAGS}>"
)

# Per-thread default stream option see https://docs.nvidia.com/cuda/cuda-runtime-api/stream-sync-behavior.html
# The per-thread default stream does not synchronize with other streams
target_compile_definitions(test_connected_components PUBLIC CUDA_API_PER_THREAD_DEFAULT_STREAM)

target_link_libraries(test_connected_components
PUBLIC
    cugraphtestutil
    mmio
    #$<BUILD_INTERFACE:conda_env>
    cugraph::cugraph
    CUDA::cublas
    CUDA::cusparse
    CUDA::cusolver
    CUDA::curand
    NCCL::NCCL
)


add_executable(test_mmio_read test_mmio_read.cpp)

set_target_properties(test_mmio_read
PROPERTIES  BUILD_RPATH    "\$ORIGIN"
            INSTALL_RPATH  "\$ORIGIN"
            CXX_STANDARD                        17
            CXX_STANDARD_REQUIRED               ON
            CUDA_STANDARD                       17
            CUDA_STANDARD_REQUIRED              ON
            INTERFACE_POSITION_INDEPENDENT_CODE ON
)

target_compile_options(test_mmio_read
    PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:${CUGRAPH_CXX_FLAGS}>"
            "$<$<COMPILE_LANGUAGE:CUDA>:${CUGRAPH_CUDA_FLAGS}>"
)

# Per-thread default stream option see https://docs.nvidia.com/cuda/cuda-runtime-api/stream-sync-behavior.html
# The per-thread default stream does not synchronize with other streams
target_compile_definitions(test_mmio_read PUBLIC CUDA_API_PER_THREAD_DEFAULT_STREAM)


target_link_libraries(test_mmio_read
PUBLIC
    cugraphtestutil
    mmio
    #$<BUILD_INTERFACE:conda_env>
    cugraph::cugraph
    CUDA::cublas
    CUDA::cusparse
    CUDA::cusolver
    CUDA::curand
    NCCL::NCCL
)
