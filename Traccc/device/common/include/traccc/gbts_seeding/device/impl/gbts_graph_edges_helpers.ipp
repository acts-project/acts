/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/utils/trigonometric_helpers.hpp"

namespace traccc::device::details {

TRACCC_HOST_DEVICE inline gbts_phi_window gbts_make_phi_window(const float lo,
                                                               const float hi) {
  constexpr float below = -traccc::device::PI_F - 1.0f;
  constexpr float above = traccc::device::PI_F + 1.0f;
  gbts_phi_window w{lo, hi, 1.0f, 0.0f, false};
  if (hi - lo >= traccc::device::TWO_PI_F) {
    w.whole = true;
  } else if (lo < -traccc::device::PI_F) {
    w = {below, hi, lo + traccc::device::TWO_PI_F, above, false};
  } else if (hi > traccc::device::PI_F) {
    w = {below, hi - traccc::device::TWO_PI_F, lo, above, false};
  }
  return w;
}

/// Binary search for the first index in [begin, end) with phis[index] >= value.
template <typename vector_t>
TRACCC_HOST_DEVICE inline unsigned int gbts_phi_lower_bound(
    const vector_t& phis, unsigned int begin, unsigned int end,
    const float value) {
  while (begin < end) {
    const unsigned int mid = begin + (end - begin) / 2u;
    if (phis[mid] < value) {
      begin = mid + 1u;
    } else {
      end = mid;
    }
  }
  return begin;
}

/// The geometry of the candidate edge node1 -> node2 as (tau, dr, dz,
/// curvature): slope dz / dr, radial and longitudinal separation and the
/// phi difference @c dphi (wrapped into [-pi, pi)) over dr. Node params
/// are (tau_min, tau_max, r, z).
TRACCC_HOST_DEVICE inline float4 gbts_make_edge_geometry(const float4 np1,
                                                         const float4 np2,
                                                         const float dphi) {
  const float dr = np2.z - np1.z;
  const float dz = np2.w - np1.w;
  return float4{dz / dr, dr, dz, dphi / dr};
}

TRACCC_HOST_DEVICE inline bool gbts_edge_passes_cuts(
    const float4 np1, const float4 np2, const float4 geo, const float dphi,
    const float delta_phi, const gbts_make_graph_edges_params& cuts) {
  if (geo.y < cuts.minDeltaRadius) {
    return false;
  }
  const float ftau = math::fabs(geo.x);
  if ((ftau < np2.x) || (ftau > np2.y) || (ftau < np1.x) || (ftau > np1.y)) {
    return false;
  }
  const float z0 = np1.w - np1.z * geo.x;
  if ((z0 < cuts.min_z0) || (z0 > cuts.max_z0)) {
    return false;
  }
  const float z_outer = z0 + cuts.maxOuterRadius * geo.x;
  if ((z_outer < cuts.cut_zMinU) || (z_outer > cuts.cut_zMaxU)) {
    return false;
  }
  if (math::fabs(dphi) > delta_phi) {
    return false;
  }
  const float curv_max = (ftau < cuts.max_Kappa_change_tau)
                             ? cuts.max_Kappa_low_tau
                             : cuts.max_Kappa_high_tau;
  if (math::fabs(geo.w) > curv_max) {
    return false;
  }
  return true;
}

/// Test the outer nodes [begin, end) with a phi in [lo, hi] against the cuts
/// of inner node (np1, phi1) and hand every accepted one to @c accept(node2,
/// np2, phi2, geometry) which is given by the caller.
template <typename accept_t>
TRACCC_HOST_DEVICE inline void gbts_walk_interval(
    const vecmem::device_vector<const float>& phis,
    const vecmem::device_vector<const float4>& packs, const unsigned int begin,
    const unsigned int end, const float lo, const float hi, const float4 np1,
    const float phi1, const float delta_phi,
    const gbts_make_graph_edges_params& cuts, accept_t& accept) {
  if ((begin >= end) || (hi < phis[begin]) || (lo > phis[end - 1u])) {
    return;
  }
  for (unsigned int j = gbts_phi_lower_bound(phis, begin, end, lo); j < end;
       ++j) {
    const float phi2 = phis[j];
    if (phi2 > hi) {
      break;
    }
    const float4 np2 = packs[j];
    const float dphi = traccc::detail::wrap_phi(phi2 - phi1);
    const float4 geo = gbts_make_edge_geometry(np1, np2, dphi);
    if (gbts_edge_passes_cuts(np1, np2, geo, dphi, delta_phi, cuts)) {
      accept(j, np2, phi2, geo);
    }
  }
}

/// Walk the delta-phi window of inner node (np1, phi1) through the outer
/// nodes [begin, end).
template <typename accept_t>
TRACCC_HOST_DEVICE inline void gbts_walk_window(
    const vecmem::device_vector<const float>& phis,
    const vecmem::device_vector<const float4>& packs, const unsigned int begin,
    const unsigned int end, const float4 np1, const float phi1,
    const float delta_phi, const gbts_make_graph_edges_params& cuts,
    accept_t& accept) {
  const float half_width = delta_phi + gbts_phi_window_eps;
  const gbts_phi_window w =
      gbts_make_phi_window(phi1 - half_width, phi1 + half_width);
  if (w.whole) {
    gbts_walk_interval(phis, packs, begin, end, -traccc::device::PI_F - 1.0f,
                       traccc::device::PI_F + 1.0f, np1, phi1, delta_phi, cuts,
                       accept);
  } else {
    // The second piece exists only when the window crosses +/-pi.
    const bool wraps = w.lo_b <= w.hi_b;
    gbts_walk_interval(phis, packs, begin, end, w.lo_a, w.hi_a, np1, phi1,
                       delta_phi, cuts, accept);
    // Wrapping gives a split window, which is handled by a "second" walk.
    if (wraps) {
      gbts_walk_interval(phis, packs, begin, end, w.lo_b, w.hi_b, np1, phi1,
                         delta_phi, cuts, accept);
    }
  }
}

}  // namespace traccc::device::details
