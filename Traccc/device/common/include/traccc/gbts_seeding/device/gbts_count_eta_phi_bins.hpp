/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Global Event Data) Payload for the @c
/// traccc::device::gbts_count_eta_phi_bins function
struct gbts_count_eta_phi_bins_payload {
  /// Number of eta bins
  unsigned int nEtaBins;
  /// Number of phi bins per eta slice
  unsigned int nPhiBins;
  /// (eta, phi) histogram of node counts
  vecmem::data::vector_view<const unsigned int> eta_phi_histo;
  /// Output: per-eta total node count (sum over phi)
  vecmem::data::vector_view<unsigned int> eta_node_counter;
  /// Output: per-eta phi prefix-sum scratch (in/out for the next kernel)
  vecmem::data::vector_view<unsigned int> phi_cusums;
};

/// @brief Sum the (eta, phi) histogram across phi for each eta bin.
///
/// One thread per eta-bin walks its phi row, writes the running per-phi sum
/// into the phi prefix-sum scratch (later turned cumulative), and stores the
/// total for that eta in the per-eta node counter.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_count_eta_phi_bins(
    const thread_id_t& thread_id,
    const gbts_count_eta_phi_bins_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_count_eta_phi_bins.ipp"
