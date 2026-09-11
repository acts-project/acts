# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Project include(s).
include(traccc-compiler-options-cpp)

# Set up the "build" of the traccc::core library.
traccc_add_library( traccc_core core TYPE SHARED
  # Common definitions.
  "include/traccc/definitions/track_parametrization.hpp"
  "include/traccc/definitions/math.hpp"
  "include/traccc/definitions/primitives.hpp"
  "include/traccc/definitions/common.hpp"
  "include/traccc/definitions/qualifiers.hpp"
  # Event data model.
  "include/traccc/edm/details/container_base.hpp"
  "include/traccc/edm/details/container_element.hpp"
  "include/traccc/edm/details/device_container.hpp"
  "include/traccc/edm/details/host_container.hpp"
  "include/traccc/edm/measurement_collection.hpp"
  "include/traccc/edm/impl/measurement_collection.ipp"
  "include/traccc/edm/measurement_helpers.hpp"
  "include/traccc/edm/impl/measurement_helpers.ipp"
  "include/traccc/edm/particle.hpp"
  "include/traccc/edm/track_parameters.hpp"
  "include/traccc/edm/container.hpp"
  "include/traccc/edm/silicon_cell_collection.hpp"
  "include/traccc/edm/impl/silicon_cell_collection.ipp"
  "include/traccc/edm/silicon_cluster_collection.hpp"
  "include/traccc/edm/spacepoint_collection.hpp"
  "include/traccc/edm/impl/spacepoint_collection.ipp"
  "include/traccc/edm/spacepoint_helpers.hpp"
  "include/traccc/edm/impl/spacepoint_helpers.ipp"
  "include/traccc/edm/seed_collection.hpp"
  "include/traccc/edm/impl/seed_collection.ipp"
  "include/traccc/edm/track_state_collection.hpp"
  "include/traccc/edm/impl/track_state_collection.ipp"
  "include/traccc/edm/track_state_helpers.hpp"
  "include/traccc/edm/impl/track_state_helpers.ipp"
  "include/traccc/edm/track_fit_outcome.hpp"
  "include/traccc/edm/track_constituent_link.hpp"
  "include/traccc/edm/track_collection.hpp"
  "include/traccc/edm/impl/track_collection.ipp"
  "include/traccc/edm/track_container.hpp"
  # Magnetic field description.
  "include/traccc/bfield/magnetic_field_types.hpp"
  "include/traccc/bfield/magnetic_field.hpp"
  "include/traccc/bfield/impl/magnetic_field.ipp"
  "include/traccc/bfield/construct_const_bfield.hpp"
  "include/traccc/bfield/impl/construct_const_bfield.ipp"
  "src/bfield/construct_const_bfield.cpp"
  # Geometry description.
  "include/traccc/geometry/detector.hpp"
  "include/traccc/geometry/module_map.hpp"
  "include/traccc/geometry/geometry.hpp"
  "include/traccc/geometry/pixel_data.hpp"
  "include/traccc/geometry/detector_design_description.hpp"
  "include/traccc/geometry/detector_conditions_description.hpp"
  # Utilities.
  "include/traccc/utils/algorithm.hpp"
  "include/traccc/utils/detray_conversion.hpp"
  "include/traccc/utils/impl/detray_conversion.ipp"
  "include/traccc/utils/type_traits.hpp"
  "include/traccc/utils/memory_resource.hpp"
  "include/traccc/utils/seed_generator.hpp"
  "include/traccc/utils/subspace.hpp"
  "include/traccc/utils/logging.hpp"
  "include/traccc/utils/prob.hpp"
  "src/utils/logging.cpp"
  # Clusterization algorithmic code.
  "include/traccc/clusterization/details/sparse_ccl.hpp"
  "include/traccc/clusterization/impl/sparse_ccl.ipp"
  "include/traccc/clusterization/sparse_ccl_algorithm.hpp"
  "src/clusterization/sparse_ccl_algorithm.cpp"
  "include/traccc/clusterization/details/measurement_creation.hpp"
  "include/traccc/clusterization/impl/measurement_creation.ipp"
  "include/traccc/clusterization/measurement_creation_algorithm.hpp"
  "src/clusterization/measurement_creation_algorithm.cpp"
  "include/traccc/clusterization/measurement_sorting_algorithm.hpp"
  "src/clusterization/measurement_sorting_algorithm.cpp"
  "include/traccc/clusterization/clusterization_algorithm.hpp"
  "src/clusterization/clusterization_algorithm.cpp"
  # Finding algorithmic code
  "include/traccc/finding/actors/ckf_aborter.hpp"
  "include/traccc/finding/actors/measurement_kalman_updater.hpp"
  "include/traccc/finding/details/combinatorial_kalman_filter_types.hpp"
  "include/traccc/finding/details/combinatorial_kalman_filter.hpp"
  "include/traccc/finding/details/progressive_kalman_filter_types.hpp"
  "include/traccc/finding/details/progressive_kalman_filter.hpp"
  "include/traccc/finding/details/run_progressive_kalman_filter.hpp"
  "include/traccc/finding/candidate_link.hpp"
  "include/traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
  "include/traccc/finding/finding_config.hpp"
  "include/traccc/finding/measurement_selector.hpp"
  "src/finding/combinatorial_kalman_filter_algorithm.cpp"
  # Fitting algorithmic code
  "include/traccc/fitting/kalman_filter/gain_matrix_updater.hpp"
  "include/traccc/fitting/kalman_filter/kalman_actor.hpp"
  "include/traccc/fitting/kalman_filter/kalman_fitter.hpp"
  "include/traccc/fitting/kalman_filter/kalman_step_aborter.hpp"
  "include/traccc/fitting/kalman_filter/statistics_updater.hpp"
  "include/traccc/fitting/kalman_filter/two_filters_smoother.hpp"
  "include/traccc/fitting/details/kalman_fitting_types.hpp"
  "include/traccc/fitting/details/kalman_fitting.hpp"
  "include/traccc/fitting/kalman_smoother.hpp"
  "include/traccc/fitting/kalman_fitting_algorithm.hpp"
  "include/traccc/fitting/triplet_fit/triplet_fitter.hpp"
  "include/traccc/fitting/triplet_fitting_algorithm.hpp"
  "src/fitting/kalman_fitting_algorithm.cpp"
  "src/fitting/triplet_fitting_algorithm.cpp"
  # Seed finding algorithmic code.
  "include/traccc/seeding/detail/lin_circle.hpp"
  "include/traccc/seeding/detail/doublet.hpp"
  "include/traccc/seeding/detail/triplet.hpp"
  "include/traccc/seeding/detail/singlet.hpp"
  "include/traccc/seeding/detail/seeding_config.hpp"
  "include/traccc/seeding/detail/spacepoint_grid.hpp"
  "include/traccc/seeding/seed_selecting_helper.hpp"
  "src/seeding/seed_filtering.hpp"
  "src/seeding/seed_filtering.cpp"
  "include/traccc/seeding/seeding_algorithm.hpp"
  "src/seeding/seeding_algorithm.cpp"
  "include/traccc/seeding/track_params_estimation_helper.hpp"
  "include/traccc/seeding/doublet_finding_helper.hpp"
  "include/traccc/seeding/spacepoint_binning_helper.hpp"
  "include/traccc/seeding/track_params_estimation.hpp"
  "src/seeding/track_params_estimation.cpp"
  "include/traccc/seeding/triplet_finding_helper.hpp"
  "src/seeding/doublet_finding.hpp"
  "src/seeding/triplet_finding.hpp"
  "include/traccc/seeding/detail/seed_finding.hpp"
  "src/seeding/seed_finding.cpp"
  "include/traccc/seeding/detail/spacepoint_binning.hpp"
  "src/seeding/spacepoint_binning.cpp"
  "include/traccc/seeding/detail/spacepoint_formation.hpp"
  "include/traccc/seeding/impl/spacepoint_formation.ipp"
  "src/seeding/silicon_pixel_spacepoint_formation.hpp"
  "include/traccc/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
  "src/seeding/silicon_pixel_spacepoint_formation_algorithm.cpp"
  #gbts seed finding config
  "include/traccc/gbts_seeding/gbts_seeding_config.hpp"
  "src/gbts_seeding/gbts_seeding_config.cpp"
  # Ambiguity resolution
  "include/traccc/ambiguity_resolution/ambiguity_resolution_config.hpp"
  "include/traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
  "src/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.cpp"
  "include/traccc/ambiguity_resolution/legacy/greedy_ambiguity_resolution_algorithm.hpp"
  "src/ambiguity_resolution/legacy/greedy_ambiguity_resolution_algorithm.cpp"
)
target_link_libraries(
    traccc_core
    PUBLIC
        vecmem::core
        covfie::core
        detray::core_array
        detray::detectors
        Acts::Core
)

string(REPLACE ";" ", " TRACCC_DETECTOR_TYPES "${TRACCC_SUPPORTED_DETECTORS}")
message(STATUS "Building with detector types: ${TRACCC_DETECTOR_TYPES}")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/include/traccc/geometry/detector_type_list.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/traccc/geometry/detector_type_list.hpp"
    @ONLY
)

target_include_directories(
    traccc_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include>
)

install(
    DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

# Any client that wants to use traccc::core with nvcc, needs to use the
# "--expt-relaxed-constexpr" flag. Otherwise nvcc has a tendency to generate
# invalid code.
target_compile_options(
    traccc_core
    INTERFACE
        $<$<AND:$<CUDA_COMPILER_ID:NVIDIA>,$<COMPILE_LANGUAGE:CUDA>>:--expt-relaxed-constexpr>
)
