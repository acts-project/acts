# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Enable CUDA as a language.
enable_language( CUDA )

# CMake 3.26 is the first to support C++20 in CUDA.
cmake_minimum_required( VERSION 3.26 )

# Find the CUDA toolkit
find_package( CUDAToolkit REQUIRED )

# Project include(s).
include( traccc-compiler-options-cpp )
include( traccc-compiler-options-cuda )

# Set up the build of the traccc::cuda library.
traccc_add_library( traccc_cuda cuda TYPE SHARED
  # Utility definitions.
  "include/traccc/cuda/utils/algorithm_base.hpp"
  "src/utils/algorithm_base.cpp"
  "src/utils/barrier.hpp"
  "src/utils/cuda_error_handling.hpp"
  "src/utils/cuda_error_handling.cpp"
  "src/utils/global_index.hpp"
  "src/utils/magnetic_field_types.hpp"
  "src/utils/sync.cuh"
  "src/utils/thread_id.hpp"
  "src/utils/utils.hpp"
  "src/utils/utils.cpp"
  "src/sanity/contiguous_on.cuh"
  "src/sanity/ordered_on.cuh"
  # Seed finding code.
  "include/traccc/cuda/seeding/seed_parameter_estimation_algorithm.hpp"
  "src/seeding/seed_parameter_estimation_algorithm.cu"
  "include/traccc/cuda/seeding/triplet_seeding_algorithm.hpp"
  "src/seeding/triplet_seeding_algorithm.cu"
  "include/traccc/cuda/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
  "src/seeding/silicon_pixel_spacepoint_formation_algorithm.cu"
  #GBTS seed finding code
  "include/traccc/cuda/gbts_seeding/gbts_seeding_algorithm.hpp"
  "src/gbts_seeding/gbts_seeding_algorithm.cu"
  "src/gbts_seeding/kernels/GbtsNodesMakingKernels.cuh"
  "src/gbts_seeding/kernels/GbtsGraphMakingKernels.cuh"
  "src/gbts_seeding/kernels/GbtsGraphProcessingKernels.cuh"
  # Clusterization
  "include/traccc/cuda/clusterization/clusterization_algorithm.hpp"
  "src/clusterization/clusterization_algorithm.cu"
  "include/traccc/cuda/clusterization/measurement_sorting_algorithm.hpp"
  "src/clusterization/measurement_sorting_algorithm.cu"
  "src/clusterization/kernels/ccl_kernel.cu"
  "src/clusterization/kernels/ccl_kernel.cuh"
  "src/clusterization/kernels/reify_cluster_data.cu"
  # Track finding algorithm(s).
  "include/traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"
  "src/finding/combinatorial_kalman_filter_algorithm.cpp"
  "src/finding/combinatorial_kalman_filter_algorithm.cu"
  "src/finding/combinatorial_kalman_filter.cuh"
  "src/finding/kernels/fill_finding_propagation_sort_keys.cu"
  "src/finding/kernels/fill_finding_propagation_sort_keys.cuh"
  "src/finding/kernels/fill_finding_duplicate_removal_sort_keys.cu"
  "src/finding/kernels/remove_duplicates.cu"
  "src/finding/kernels/build_tracks.cu"
  "src/finding/kernels/build_tracks.cuh"
  "src/finding/kernels/find_tracks.cuh"
  "src/finding/kernels/gather_best_tips_per_measurement.cu"
  "src/finding/kernels/gather_measurement_votes.cu"
  "src/finding/kernels/update_tip_length_buffer.cu"
  "src/finding/kernels/propagate_to_next_surface.hpp"
  # Ambiguity resolution
  "include/traccc/cuda/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
  "src/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.cu"
  "src/ambiguity_resolution/kernels/add_block_offset.cu"
  "src/ambiguity_resolution/kernels/add_block_offset.cuh"
  "src/ambiguity_resolution/kernels/block_inclusive_scan.cu"
  "src/ambiguity_resolution/kernels/block_inclusive_scan.cuh"
  "src/ambiguity_resolution/kernels/count_shared_measurements.cu"
  "src/ambiguity_resolution/kernels/count_shared_measurements.cuh"
  "src/ambiguity_resolution/kernels/fill_inverted_ids.cu"
  "src/ambiguity_resolution/kernels/fill_inverted_ids.cuh"
  "src/ambiguity_resolution/kernels/fill_tracks_per_measurement.cu"
  "src/ambiguity_resolution/kernels/fill_tracks_per_measurement.cuh"
  "src/ambiguity_resolution/kernels/fill_unique_meas_id_map.cu"
  "src/ambiguity_resolution/kernels/fill_unique_meas_id_map.cuh"
  "src/ambiguity_resolution/kernels/fill_track_candidates.cu"
  "src/ambiguity_resolution/kernels/fill_track_candidates.cuh"
  "src/ambiguity_resolution/kernels/fill_vectors.cu"
  "src/ambiguity_resolution/kernels/fill_vectors.cuh"
  "src/ambiguity_resolution/kernels/update_status.cu"
  "src/ambiguity_resolution/kernels/update_status.cuh"
  "src/ambiguity_resolution/kernels/rearrange_tracks.cu"
  "src/ambiguity_resolution/kernels/rearrange_tracks.cuh"
  "src/ambiguity_resolution/kernels/scan_block_offsets.cu"
  "src/ambiguity_resolution/kernels/scan_block_offsets.cuh"
  "src/ambiguity_resolution/kernels/sort_tracks_per_measurement.cu"
  "src/ambiguity_resolution/kernels/sort_tracks_per_measurement.cuh"
  "src/ambiguity_resolution/kernels/sort_updated_tracks.cu"
  "src/ambiguity_resolution/kernels/sort_updated_tracks.cuh"
  "src/ambiguity_resolution/kernels/remove_tracks.cu"
  "src/ambiguity_resolution/kernels/remove_tracks.cuh"
  # Track fitting algorithm(s).
  "include/traccc/cuda/fitting/kalman_fitting_algorithm.hpp"
  "src/fitting/kalman_fitting_algorithm.cpp"
  "src/fitting/kalman_fitting_algorithm.cu"
  "src/fitting/kalman_fitting.cuh"
  "src/fitting/kernels/fill_fitting_sort_keys.cu"
  "src/fitting/kernels/fit_prelude.cu"
)

set(TRACCC_CUDA_SUPPORTED_BFIELDS "const;inhom_global;inhom_texture")

set(KERNEL_SPECIALIZATION_PY "${PROJECT_SOURCE_DIR}/codegen/kernel_specialization/gen_kernel_specialization.py")

# Generate specializations of find_tracks
foreach(DETECTOR_NAME ${TRACCC_SUPPORTED_DETECTORS})
    set(GENERATED_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/src/finding/kernels/specializations/find_tracks_${DETECTOR_NAME}.cu")
    set(TEMPLATE_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/finding/kernels/specializations/find_tracks.cu.template")
    add_custom_command(
        OUTPUT "${GENERATED_SOURCE}"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/src/finding/kernels/specializations/"
        COMMAND
            Python::Interpreter
            "${KERNEL_SPECIALIZATION_PY}"
            "${TEMPLATE_SOURCE}"
            -o "${GENERATED_SOURCE}"
            --detector "${DETECTOR_NAME}"
        DEPENDS "${KERNEL_SPECIALIZATION_PY}" "${TEMPLATE_SOURCE}"
        COMMENT "Generating kernel specialization of `find_tracks` for ${DETECTOR_NAME}"
    )
    target_sources(traccc_cuda PRIVATE ${GENERATED_SOURCE})
endforeach()

# Generate specializations for propagate_to_next_surface
foreach(DETECTOR_NAME ${TRACCC_SUPPORTED_DETECTORS})
    foreach(BFIELD_NAME ${TRACCC_CUDA_SUPPORTED_BFIELDS})
        set(GENERATED_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/src/finding/kernels/specializations/propagate_to_next_surface_${DETECTOR_NAME}_${BFIELD_NAME}.cu")
        set(TEMPLATE_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/finding/kernels/specializations/propagate_to_next_surface.cu.template")
        add_custom_command(
            OUTPUT "${GENERATED_SOURCE}"
            COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/src/finding/kernels/specializations/"
            COMMAND
                Python::Interpreter
                "${KERNEL_SPECIALIZATION_PY}"
                "${TEMPLATE_SOURCE}"
                -o "${GENERATED_SOURCE}"
                --detector "${DETECTOR_NAME}"
                --bfield "${BFIELD_NAME}"
                --model cuda
            DEPENDS "${KERNEL_SPECIALIZATION_PY}" "${TEMPLATE_SOURCE}"
            COMMENT "Generating kernel specialization of `propagate_to_next_surface` for ${DETECTOR_NAME} and ${BFIELD_NAME}"
        )
        target_sources(traccc_cuda PRIVATE ${GENERATED_SOURCE})
    endforeach()
endforeach()

# Generate specializations for fit_forward
foreach(DETECTOR_NAME ${TRACCC_SUPPORTED_DETECTORS})
    foreach(BFIELD_NAME ${TRACCC_CUDA_SUPPORTED_BFIELDS})
        set(GENERATED_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/src/fitting/kernels/specializations/fit_forward_${DETECTOR_NAME}_${BFIELD_NAME}.cu")
        set(TEMPLATE_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/fitting/kernels/specializations/fit_forward.cu.template")
        add_custom_command(
            OUTPUT "${GENERATED_SOURCE}"
            COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/src/fitting/kernels/specializations/"
            COMMAND
                Python::Interpreter
                "${KERNEL_SPECIALIZATION_PY}"
                "${TEMPLATE_SOURCE}"
                -o "${GENERATED_SOURCE}"
                --detector "${DETECTOR_NAME}"
                --bfield "${BFIELD_NAME}"
                --model cuda
            DEPENDS "${KERNEL_SPECIALIZATION_PY}" "${TEMPLATE_SOURCE}"
            COMMENT "Generating kernel specialization of `fit_forward` for ${DETECTOR_NAME} and ${BFIELD_NAME}"
        )
        target_sources(traccc_cuda PRIVATE ${GENERATED_SOURCE})
    endforeach()
endforeach()

# Generate specializations for fit_backward
foreach(DETECTOR_NAME ${TRACCC_SUPPORTED_DETECTORS})
    foreach(BFIELD_NAME ${TRACCC_CUDA_SUPPORTED_BFIELDS})
        set(GENERATED_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/src/fitting/kernels/specializations/fit_backward_${DETECTOR_NAME}_${BFIELD_NAME}.cu")
        set(TEMPLATE_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/fitting/kernels/specializations/fit_backward.cu.template")
        add_custom_command(
            OUTPUT "${GENERATED_SOURCE}"
            COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/src/fitting/kernels/specializations/"
            COMMAND
                Python::Interpreter
                "${KERNEL_SPECIALIZATION_PY}"
                "${TEMPLATE_SOURCE}"
                -o "${GENERATED_SOURCE}"
                --detector "${DETECTOR_NAME}"
                --bfield "${BFIELD_NAME}"
                --model cuda
            DEPENDS "${KERNEL_SPECIALIZATION_PY}" "${TEMPLATE_SOURCE}"
            COMMENT "Generating kernel specialization of `fit_backward` for ${DETECTOR_NAME} and ${BFIELD_NAME}"
        )
        target_sources(traccc_cuda PRIVATE ${GENERATED_SOURCE})
    endforeach()
endforeach()

if(TRACCC_ENABLE_NVTX_PROFILING)
    traccc_add_library(
        nvtx_instrumentation
        TYPE OBJECT
        "src/inst_nvtx.cpp"
    )

    set_target_properties(
        nvtx_instrumentation
        PROPERTIES
        POSITION_INDEPENDENT_CODE ON
    )

    target_link_libraries(
        nvtx_instrumentation
        PRIVATE CUDA::nvToolsExt ${CMAKE_DL_LIBS}
    )

    target_compile_options(
        traccc_cuda
        PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-finstrument-functions>
        $<$<COMPILE_LANGUAGE:CXX>:-finstrument-functions>)

    target_link_libraries(
        traccc_cuda
        PRIVATE
        nvtx_instrumentation
    )
endif()

target_link_libraries( traccc_cuda
  PUBLIC traccc::core traccc::cuda_utils detray::core vecmem::core
  PRIVATE CUDA::cudart traccc::device_common vecmem::cuda covfie::cuda )

# Set up Thrust specifically for the traccc::cuda library.
thrust_create_target( traccc::cuda_thrust
  HOST CPP
  DEVICE CUDA )
target_link_libraries( traccc::cuda_thrust INTERFACE CUDA::cudart )
target_link_libraries( traccc_cuda PRIVATE traccc::cuda_thrust )

# For CUDA 11 turn on separable compilation. This is necessary for using
# Thrust 2.1.0.
if( ( "${CMAKE_CUDA_COMPILER_ID}" STREQUAL "NVIDIA" ) AND
    ( "${CMAKE_CUDA_COMPILER_VERSION}" VERSION_LESS "12.0" ) )
  set_target_properties( traccc_cuda PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON )
endif()
