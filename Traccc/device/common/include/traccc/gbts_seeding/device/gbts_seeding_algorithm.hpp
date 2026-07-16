/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/algorithm_base.hpp"
#include "traccc/gbts_seeding/device/gbts_add_terminus_to_path_store.hpp"
#include "traccc/gbts_seeding/device/gbts_bid_seeds_for_hits.hpp"
#include "traccc/gbts_seeding/device/gbts_bin_spacepoints.hpp"
#include "traccc/gbts_seeding/device/gbts_compress_graph.hpp"
#include "traccc/gbts_seeding/device/gbts_convert_seeds.hpp"
#include "traccc/gbts_seeding/device/gbts_count_eta_phi_bins.hpp"
#include "traccc/gbts_seeding/device/gbts_count_spacepoints_by_layer.hpp"
#include "traccc/gbts_seeding/device/gbts_count_terminus_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_fill_path_store.hpp"
#include "traccc/gbts_seeding/device/gbts_find_minmax_radius.hpp"
#include "traccc/gbts_seeding/device/gbts_fit_segments.hpp"
#include "traccc/gbts_seeding/device/gbts_link_graph_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_make_graph_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_match_graph_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_prefix_sum_eta_phi_bins.hpp"
#include "traccc/gbts_seeding/device/gbts_rebid_seeds_for_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_reindex_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_reset_edge_bids.hpp"
#include "traccc/gbts_seeding/device/gbts_run_cca_iteration.hpp"
#include "traccc/gbts_seeding/device/gbts_sort_nodes.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/vector.hpp>

// System include(s).
#include <cstdint>
#include <memory>
#include <utility>

namespace traccc::device {

/// @brief Main algorithm for performing GBTS seeding on a device
/// (backend-agnostic).
///
/// The algorithm orchestrates the sequence of kernel launches and host-side
/// synchronisations.  Backend-specific subclasses are responsible for
/// implementing the individual kernel launchers.
///
/// This algorithm returns a buffer which is not necessarily filled yet. A
/// synchronisation statement is required before destroying this buffer.
///
class gbts_seeding_algorithm
    : public algorithm<edm::seed_collection::buffer(
          const edm::spacepoint_collection::const_view&,
          const edm::measurement_collection::const_view&)>,
      public messaging,
      public algorithm_base {
 public:
  /// Constructor for the GBTS seed finding algorithm
  ///
  /// @param cfg The GBTS seed finding configuration
  /// @param mr The memory resource(s) to use in the algorithm
  /// @param copy The copy object to use for copying data between device
  ///             and host memory blocks
  /// @param logger The logger instance to use
  ///
  gbts_seeding_algorithm(
      const gbts_seedfinder_config& cfg, const memory_resource& mr,
      const vecmem::copy& copy,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone());

  /// Destructor
  virtual ~gbts_seeding_algorithm() = default;

  /// Operator executing the algorithm.
  ///
  /// @param spacepoints is a view of all spacepoints in the event
  /// @param measurements is a view of all measurements in the event
  /// @return the buffer of track seeds reconstructed from the spacepoints
  ///
  output_type operator()(
      const edm::spacepoint_collection::const_view& spacepoints,
      const edm::measurement_collection::const_view& measurements)
      const override;

 protected:
  /// @name Kernel launchers (to be implemented by backends)
  ///
  /// Each launcher receives the payload of the device function it runs;
  /// the payload types are defined in the per-function device headers.
  /// @{

  /// Spacepoint-by-layer counting kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_count_spacepoints_by_layer_kernel(
      const gbts_count_spacepoints_by_layer_payload& payload) const = 0;

  /// Spacepoint-binning kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_bin_spacepoints_kernel(
      const gbts_bin_spacepoints_payload& payload) const = 0;

  /// Eta-phi counting kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_count_eta_phi_bins_kernel(
      const gbts_count_eta_phi_bins_payload& payload) const = 0;

  /// Eta-phi prefix-sum kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_prefix_sum_eta_phi_bins_kernel(
      const gbts_prefix_sum_eta_phi_bins_payload& payload) const = 0;

  /// Node sorting kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_sort_nodes_kernel(
      const gbts_sort_nodes_payload& payload) const = 0;

  /// Min/max radius per eta-bin kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_find_minmax_radius_kernel(
      const gbts_find_minmax_radius_payload& payload) const = 0;

  /// Graph edge-making kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_make_graph_edges_kernel(
      const gbts_make_graph_edges_payload& payload) const = 0;

  /// Graph edge-linking kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_link_graph_edges_kernel(
      const gbts_link_graph_edges_payload& payload) const = 0;

  /// Graph edge-matching kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_match_graph_edges_kernel(
      const gbts_match_graph_edges_payload& payload) const = 0;

  /// Edge re-indexing kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_reindex_edges_kernel(
      const gbts_reindex_edges_payload& payload) const = 0;

  /// Graph compression kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_compress_graph_kernel(
      const gbts_compress_graph_payload& payload) const = 0;

  /// CCA (connected-components iteration) kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_run_cca_iteration_kernel(
      const gbts_run_cca_iteration_payload& payload) const = 0;

  /// Terminus-edge counting kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_count_terminus_edges_kernel(
      const gbts_count_terminus_edges_payload& payload) const = 0;

  /// Terminus-to-path-store seeding kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_add_terminus_to_path_store_kernel(
      const gbts_add_terminus_to_path_store_payload& payload) const = 0;

  /// Path-store-filling kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_fill_path_store_kernel(
      const gbts_fill_path_store_payload& payload) const = 0;

  /// Segment fitting kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_fit_segments_kernel(
      const gbts_fit_segments_payload& payload) const = 0;

  /// Edge-bid reset kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_reset_edge_bids_kernel(
      const gbts_reset_edge_bids_payload& payload) const = 0;

  /// Edge re-bid kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_rebid_seeds_for_edges_kernel(
      const gbts_rebid_seeds_for_edges_payload& payload) const = 0;

  /// Seeds-bid-for-hits kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_bid_seeds_for_hits_kernel(
      const gbts_bid_seeds_for_hits_payload& payload) const = 0;

  /// GBTS seed conversion kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void gbts_convert_seeds_kernel(
      const gbts_convert_seeds_payload& payload) const = 0;

  /// @}

 private:
  /// @name Pipeline stages
  ///
  /// The pipeline is split into three stage methods so that each stage's
  /// transient device buffers are local and freed as soon as the stage
  /// returns; only the cross-stage handles below survive between stages.
  /// @{

  /// Outputs of the node-making stage that are consumed downstream.
  struct node_making_output {
    /// Reduced (x, y, z, w) per original spacepoint (used by seed
    /// extraction)
    vecmem::data::vector_buffer<float4> reducedSP;
    /// Per-node (tau_min, tau_max, r, z) (used by graph making)
    vecmem::data::vector_buffer<float4> node_params;
    /// Per-node phi (used by graph making)
    vecmem::data::vector_buffer<float> node_phi;
    /// Per-sorted-slot original spacepoint index (used by graph making)
    vecmem::data::vector_buffer<unsigned int> node_index;
    /// Per-eta (rmin, rmax) pair, host (used by graph making)
    vecmem::vector<float> bin_rads;
    /// Per-eta (begin, end) node ranges, host (used by graph making)
    vecmem::vector<unsigned int> eta_bin_views;
    /// Number of GBTS nodes (0 == nothing to do)
    unsigned int nNodes = 0;
  };

  /// Outputs of the graph-making stage that are consumed by seed extraction.
  struct graph_making_output {
    /// Compacted, row-major graph
    vecmem::data::vector_buffer<unsigned int> output_graph;
    /// Number of edges that survived re-indexing (0 == nothing to do)
    unsigned int nConnectedEdges = 0;
  };

  /// Stage 1: count, bin, sort and characterise nodes.
  node_making_output make_nodes(
      const edm::spacepoint_collection::const_view& spacepoints,
      const edm::measurement_collection::const_view& measurements) const;

  /// Stage 2: build, link, match and compress the edge graph. The per-node
  /// buffers are taken by value so they are released when this stage returns.
  graph_making_output create_edges(
      vecmem::data::vector_buffer<float4> node_params,
      vecmem::data::vector_buffer<float> node_phi,
      vecmem::data::vector_buffer<unsigned int> node_index,
      const vecmem::vector<float>& bin_rads,
      const vecmem::vector<unsigned int>& eta_bin_views,
      const unsigned int nNodes,
      vecmem::data::vector_buffer<unsigned int>& counters_buf,
      vecmem::vector<unsigned int>& h_counters) const;

  /// Stage 3: run the CCA, extract paths, fit and disambiguate into seeds.
  edm::seed_collection::buffer extract_seeds(
      vecmem::data::vector_buffer<unsigned int>& output_graph,
      vecmem::data::vector_buffer<float4>& reducedSP,
      const unsigned int nConnectedEdges, const unsigned int nSp,
      vecmem::data::vector_buffer<unsigned int>& counters_buf,
      vecmem::vector<unsigned int>& h_counters) const;

  /// @}

  /// GBTS seed-finding configuration.
  gbts_seedfinder_config m_config;

};  // class gbts_seeding_algorithm

}  // namespace traccc::device
