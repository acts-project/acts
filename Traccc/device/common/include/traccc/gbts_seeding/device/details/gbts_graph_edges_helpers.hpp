/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::device::details {

/// Padding of the phi search window; the exact delta-phi cut decides.
inline constexpr float gbts_phi_window_eps = 1e-5f;

/// A phi interval [lo, hi] split at +/-pi into up to two pieces, both
/// ascending in phi. The second piece is unused when hi_b < lo_b.
struct gbts_phi_window {
  float lo_a, hi_a;
  float lo_b, hi_b;
  bool whole;  ///< the interval covers the full circle
};

/// Split the phi interval [lo, hi] at +/-pi
TRACCC_HOST_DEVICE inline gbts_phi_window gbts_make_phi_window(const float lo,
                                                               const float hi);

/// First index in [begin, end) of the ascending @c phis with a value >= @c
/// value
template <typename vector_t>
TRACCC_HOST_DEVICE inline unsigned int gbts_phi_lower_bound(
    const vector_t& phis, unsigned int begin, unsigned int end,
    const float value);

/// The geometry of the candidate edge node1 -> node2 as (tau, dr, dz,
/// curvature).
TRACCC_HOST_DEVICE inline float4 gbts_make_edge_geometry(const float4 np1,
                                                         const float4 np2,
                                                         const float dphi);

/// The doublet cuts for the candidate edge node1 -> node2 with geometry
/// @c geo (see gbts_make_edge_geometry). Node params are (tau_min,
/// tau_max, r, z).
TRACCC_HOST_DEVICE inline bool gbts_edge_passes_cuts(
    const float4 np1, const float4 np2, const float4 geo, const float dphi,
    const float delta_phi, const gbts_make_graph_edges_params& cuts);

/// Test the outer nodes [begin, end) with a phi in [lo, hi] against the cuts
/// of inner node (np1, phi1) and hand every accepted one to @c accept(node2,
/// np2, phi2, geometry), in ascending node order.
template <typename accept_t>
TRACCC_HOST_DEVICE inline void gbts_walk_interval(
    const vecmem::device_vector<const float>& phis,
    const vecmem::device_vector<const float4>& packs, const unsigned int begin,
    const unsigned int end, const float lo, const float hi, const float4 np1,
    const float phi1, const float delta_phi,
    const gbts_make_graph_edges_params& cuts, accept_t& accept);

/// Walk the delta-phi window of inner node (np1, phi1) through the outer
/// nodes [begin, end).
template <typename accept_t>
TRACCC_HOST_DEVICE inline void gbts_walk_window(
    const vecmem::device_vector<const float>& phis,
    const vecmem::device_vector<const float4>& packs, const unsigned int begin,
    const unsigned int end, const float4 np1, const float phi1,
    const float delta_phi, const gbts_make_graph_edges_params& cuts,
    accept_t& accept);

}  // namespace traccc::device::details

#include "traccc/gbts_seeding/device/impl/gbts_graph_edges_helpers.ipp"
