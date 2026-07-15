/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/sycl/utils/algorithm_base.hpp"

// Project include(s).
#include "traccc/gbts_seeding/device/gbts_seeding_algorithm.hpp"

namespace traccc::sycl {

/// @brief Main algorithm for performing GBTS seeding using oneAPI/SYCL.
///
/// This algorithm returns a buffer which is not necessarily filled yet. A
/// synchronisation statement is required before destroying this buffer.
///
class gbts_seeding_algorithm : public device::gbts_seeding_algorithm,
                               public sycl::algorithm_base {

    public:
    /// Constructor for the GBTS seed finding algorithm
    ///
    /// @param cfg The GBTS seed finding configuration
    /// @param mr The memory resource(s) to use in the algorithm
    /// @param copy The copy object to use for copying data between device
    ///             and host memory blocks
    /// @param queue The SYCL queue to perform the operations in
    /// @param logger The logger instance to use
    ///
    gbts_seeding_algorithm(
        const gbts_seedfinder_config& cfg, const memory_resource& mr,
        const vecmem::copy& copy, queue_wrapper& queue,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    private:
    /// @name Function(s) inherited from @c
    /// traccc::device::gbts_seeding_algorithm
    /// @{

    void gbts_count_spacepoints_by_layer_kernel(
        const device::gbts_count_spacepoints_by_layer_payload& payload)
        const override;
    void gbts_bin_spacepoints_kernel(
        const device::gbts_bin_spacepoints_payload& payload) const override;
    void gbts_count_eta_phi_bins_kernel(
        const device::gbts_count_eta_phi_bins_payload& payload) const override;
    void gbts_prefix_sum_eta_phi_bins_kernel(
        const device::gbts_prefix_sum_eta_phi_bins_payload& payload)
        const override;
    void gbts_sort_nodes_kernel(
        const device::gbts_sort_nodes_payload& payload) const override;
    void gbts_find_minmax_radius_kernel(
        const device::gbts_find_minmax_radius_payload& payload) const override;
    void gbts_make_graph_edges_kernel(
        const device::gbts_make_graph_edges_payload& payload) const override;
    void gbts_link_graph_edges_kernel(
        const device::gbts_link_graph_edges_payload& payload) const override;
    void gbts_match_graph_edges_kernel(
        const device::gbts_match_graph_edges_payload& payload) const override;
    void gbts_reindex_edges_kernel(
        const device::gbts_reindex_edges_payload& payload) const override;
    void gbts_compress_graph_kernel(
        const device::gbts_compress_graph_payload& payload) const override;
    void gbts_run_cca_iteration_kernel(
        const device::gbts_run_cca_iteration_payload& payload) const override;
    void gbts_count_terminus_edges_kernel(
        const device::gbts_count_terminus_edges_payload& payload)
        const override;
    void gbts_add_terminus_to_path_store_kernel(
        const device::gbts_add_terminus_to_path_store_payload& payload)
        const override;
    void gbts_fill_path_store_kernel(
        const device::gbts_fill_path_store_payload& payload) const override;
    void gbts_fit_segments_kernel(
        const device::gbts_fit_segments_payload& payload) const override;
    void gbts_reset_edge_bids_kernel(
        const device::gbts_reset_edge_bids_payload& payload) const override;
    void gbts_rebid_seeds_for_edges_kernel(
        const device::gbts_rebid_seeds_for_edges_payload& payload)
        const override;
    void gbts_bid_seeds_for_hits_kernel(
        const device::gbts_bid_seeds_for_hits_payload& payload) const override;
    void gbts_convert_seeds_kernel(
        const device::gbts_convert_seeds_payload& payload) const override;

    /// @}

};  // class gbts_seeding_algorithm

}  // namespace traccc::sycl
