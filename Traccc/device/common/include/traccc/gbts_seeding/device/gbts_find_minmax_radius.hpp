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
/// traccc::device::gbts_find_minmax_radius function
struct gbts_find_minmax_radius_payload {
  /// Number of eta bins
  unsigned int nEtaBins;
  /// Per-eta (begin, end) node range, as 2*nEtaBins flat ints
  vecmem::data::vector_view<const unsigned int> eta_bin_views;
  /// Per-node (tau_min, tau_max, r, z) (only r is read here)
  vecmem::data::vector_view<const float4> node_params;
  /// Output: per-eta (rmin, rmax) pair, flat (2*nEtaBins floats)
  vecmem::data::vector_view<float> bin_rads;
};

/// @brief Compute the per-eta-bin minimum and maximum radius.
///
/// One thread per eta-bin scans the bin's node range and writes the (rmin,
/// rmax) pair into the output; the host uses these to estimate the maximum
/// delta-R for each bin pair.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_find_minmax_radius(
    const thread_id_t& thread_id,
    const gbts_find_minmax_radius_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_find_minmax_radius.ipp"
