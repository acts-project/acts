# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2022-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Project include(s).
include( traccc-compiler-options-cpp )

# Declare the traccc::device_common library.
traccc_add_library( traccc_device_common device_common
   TYPE SHARED
   # General code.
   "include/traccc/device/algorithm_base.hpp"
   "src/device/algorithm_base.cpp"
   "include/traccc/device/prefix_sum_element.hpp"
   "include/traccc/device/container_h2d_copy_alg.hpp"
   "include/traccc/device/impl/container_h2d_copy_alg.ipp"
   "include/traccc/device/container_d2h_copy_alg.hpp"
   "include/traccc/device/impl/container_d2h_copy_alg.ipp"
   # EDM class(es).
   "include/traccc/edm/device/seeding_global_counter.hpp"
   "include/traccc/edm/device/doublet_counter.hpp"
   "include/traccc/edm/device/triplet_counter.hpp"
   "include/traccc/edm/device/device_doublet.hpp"
   "include/traccc/edm/device/device_triplet.hpp"
   "include/traccc/edm/device/sort_key.hpp"
   "include/traccc/edm/device/identity_projector.hpp"
   # Clusterization code.
   "include/traccc/clusterization/device/reduce_problem_cell.hpp"
   "include/traccc/clusterization/device/impl/reduce_problem_cell.ipp"
   "include/traccc/clusterization/device/aggregate_cluster.hpp"
   "include/traccc/clusterization/device/impl/aggregate_cluster.ipp"
   "include/traccc/clusterization/device/geo_id_based_sorter.hpp"
   "include/traccc/clusterization/device/ccl_kernel_definitions.hpp"
   "include/traccc/clusterization/device/ccl_kernel.hpp"
   "include/traccc/clusterization/device/sorting_index_filler.hpp"
   "include/traccc/clusterization/device/clusterization_algorithm.hpp"
   "src/clusterization/clusterization_algorithm.cpp"
   # Spacepoint binning function(s).
   "include/traccc/seeding/device/count_grid_capacities.hpp"
   "include/traccc/seeding/device/impl/count_grid_capacities.ipp"
   "include/traccc/seeding/device/populate_grid.hpp"
   "include/traccc/seeding/device/impl/populate_grid.ipp"
   # Seed finding code.
   "include/traccc/seeding/device/form_spacepoints.hpp"
   "include/traccc/seeding/device/impl/form_spacepoints.ipp"
   "include/traccc/seeding/device/count_doublets.hpp"
   "include/traccc/seeding/device/impl/count_doublets.ipp"
   "include/traccc/seeding/device/find_doublets.hpp"
   "include/traccc/seeding/device/impl/find_doublets.ipp"
   "include/traccc/seeding/device/count_triplets.hpp"
   "include/traccc/seeding/device/impl/count_triplets.ipp"
   "include/traccc/seeding/device/reduce_triplet_counts.hpp"
   "include/traccc/seeding/device/impl/reduce_triplet_counts.ipp"
   "include/traccc/seeding/device/find_triplets.hpp"
   "include/traccc/seeding/device/impl/find_triplets.ipp"
   "include/traccc/seeding/device/update_triplet_weights.hpp"
   "include/traccc/seeding/device/impl/update_triplet_weights.ipp"
   "include/traccc/seeding/device/select_seeds.hpp"
   "include/traccc/seeding/device/impl/select_seeds.ipp"
   "include/traccc/seeding/device/triplet_seeding_algorithm.hpp"
   "src/seeding/triplet_seeding_algorithm.cpp"
   "include/traccc/seeding/device/silicon_pixel_spacepoint_formation_algorithm.hpp"
   "src/seeding/silicon_pixel_spacepoint_formation_algorithm.cpp"
   "include/traccc/seeding/device/seed_parameter_estimation_algorithm.hpp"
   "src/seeding/seed_parameter_estimation_algorithm.cpp"
   # Track parameters estimation function(s).
   "include/traccc/seeding/device/estimate_track_params.hpp"
   "include/traccc/seeding/device/impl/estimate_track_params.ipp"
   # Track finding funtions(s).
   "include/traccc/finding/device/geo_id_surface_comparator.hpp"
   "include/traccc/finding/device/build_tracks.hpp"
   "include/traccc/finding/device/find_tracks.hpp"
   "include/traccc/finding/device/gather_best_tips_per_measurement.hpp"
   "include/traccc/finding/device/gather_measurement_votes.hpp"
   "include/traccc/finding/device/fill_finding_propagation_sort_keys.hpp"
   "include/traccc/finding/device/propagate_to_next_surface.hpp"
   "include/traccc/finding/device/update_tip_length_buffer.hpp"
   "include/traccc/finding/device/impl/build_tracks.ipp"
   "include/traccc/finding/device/impl/find_tracks.ipp"
   "include/traccc/finding/device/impl/gather_best_tips_per_measurement.ipp"
   "include/traccc/finding/device/impl/gather_measurement_votes.ipp"
   "include/traccc/finding/device/impl/fill_finding_propagation_sort_keys.ipp"
   "include/traccc/finding/device/impl/propagate_to_next_surface.ipp"
   "include/traccc/finding/device/impl/update_tip_length_buffer.ipp"
   # Track fitting funtions(s).
   "include/traccc/fitting/device/fit.hpp"
   "include/traccc/fitting/device/fill_fitting_sort_keys.hpp"
   "include/traccc/fitting/device/impl/fill_fitting_sort_keys.ipp"
   )
target_link_libraries( traccc_device_common
   PUBLIC traccc::core vecmem::core )
