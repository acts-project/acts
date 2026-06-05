# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Declare the cpu algorithm test(s).
traccc_add_test(cpu
    "compare_with_acts_seeding.cpp"
    "seq_single_module.cpp"
    "test_ambiguity_resolution.cpp"
    "test_axis.cpp"
    "test_cca.cpp"
    "test_ckf_combinatorics_telescope.cpp"
    "test_ckf_sparse_tracks_telescope.cpp"
    "test_clusterization_resolution.cpp"
    "test_copy.cpp"
    "test_grid2.cpp"
    "test_kalman_fitter_hole_count.cpp"
    "test_kalman_fitter_momentum_resolution.cpp"
    "test_kalman_fitter_telescope.cpp"
    "test_triplet_fitter_telescope.cpp"
    "test_kalman_fitter_wire_chamber.cpp"
    "test_populator.cpp"
    "test_ranges.cpp"
    "test_seeding.cpp"
    "test_serializer.cpp"
    "test_simulation.cpp"
    "test_spacepoint_formation.cpp"
    "test_track_params_estimation.cpp"
    "test_sanity_ordered_on.cpp"
    "test_sanity_contiguous_on.cpp"
    LINK_LIBRARIES GTest::gtest_main vecmem::core
    traccc_tests_common traccc::core traccc::io traccc::performance
    traccc::simulation detray::core detray::io detray::test_common  covfie::core )
