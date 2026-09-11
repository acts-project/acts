/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/particle.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/io/csv/hit.hpp"
#include "traccc/io/csv/particle.hpp"

// Detray include(s).
#include <detray/navigation/intersection/intersection.hpp>

// Detray test include(s)
#include <detray/test/utils/inspectors.hpp>  //< candidate_record type

// System include(s)
#include <algorithm>

namespace traccc::propagation_validator {

template <typename detector_t>
using candidate_type = detray::intersection_record<detector_t>;

/// Transcribe the hits for a particle to a candidate trace for the detray
/// propagation validation tools
///
/// @param ctx the geometric context
/// @param det the detector
/// @param ptc the truth particle
/// @param hits all hits in the event
/// @param n_hits_for_particle expected number of hits for the particle
///
/// @returns a vector of intersection candidate records
template <typename detector_t>
auto transcribe_to_trace(const typename detector_t::geometry_context ctx,
                         const detector_t& det,
                         const traccc::io::csv::particle& ptc,
                         const std::vector<traccc::io::csv::hit>& hits,
                         const std::size_t n_hits_for_particle = 10u) {
  using intersection_t = typename candidate_type<detector_t>::intersection_type;

  detray::dvector<candidate_type<detector_t>> candidates{};
  candidates.reserve(n_hits_for_particle);

  // Fill the hits into the candidate trace
  for (const auto& h : hits) {
    if (h.particle_id != ptc.particle_id) {
      continue;
    }

    // Rough estimate of path
    const point3 pos{h.tx, h.ty, h.tz};
    const vector3 mom{h.tpx, h.tpy, h.tpz};
    const vector3 dir = vector::normalize(mom);
    const scalar path{vector::norm(pos)};

    // Corresponding surface
    const detray::geometry::identifier geo_id{h.geometry_id};
    const auto sf_desc = det.surfaces().at(geo_id.index());
    const auto sf = detray::tracking_surface{det, geo_id};

    // Build an intersection from the hit
    using nav_link_t = typename intersection_t::nav_link_t;
    auto loc_pos = sf.global_to_local(ctx, pos, dir);
    intersection_t intr{sf_desc,
                        path,
                        static_cast<nav_link_t>(geo_id.volume()),
                        detray::intersection::status::e_inside,
                        true,
                        loc_pos};

    candidates.emplace_back(pos, dir, intr, ptc.q, vector::norm(mom));
  }

  // Sort records by intersection distance to origin of the trajectory
  auto sort_by_path = [&](const candidate_type<detector_t>& a,
                          const candidate_type<detector_t>& b) -> bool {
    return (a.intersection < b.intersection);
  };

  std::ranges::stable_sort(candidates, sort_by_path);

  return candidates;
}

/// Transcribe the measurements for a particle to a candidate trace for the
/// detray propagation validation tools
///
/// @param ctx the geometric context
/// @param det the detector
/// @param ptc the truth particle
/// @param ptc_to_meas_map map the truth particle to its measurements
/// @param n_meas_for_particle expected number of ,easurements for the particle
///
/// @returns a vector of intersection candidate records
template <typename detector_t>
auto transcribe_to_trace(
    const typename detector_t::geometry_context ctx, const detector_t& det,
    const traccc::particle& ptc,
    const std::map<traccc::particle,
                   std::vector<edm::measurement_collection::host::object_type>>&
        ptc_to_meas_map,
    const std::size_t n_meas_for_particle = 10u) {
  using intersection_t = typename candidate_type<detector_t>::intersection_type;

  vecmem::vector<candidate_type<detector_t>> candidates{};
  candidates.reserve(n_meas_for_particle);

  // TODO: Not accurate for every measurement
  const scalar p{vector::norm(ptc.momentum)};
  const scalar q{ptc.charge};

  // Fill the hits into the candidate trace
  for (const auto& meas : ptc_to_meas_map.at(ptc)) {
    // Corresponding surface
    const detray::geometry::identifier geo_id{meas.surface_link()};
    const auto sf_desc = det.surfaces().at(geo_id.index());
    const auto sf = detray::tracking_surface{det, geo_id};

    // TODO: Use correct track direction at measurement for line sf.
    const vector3 dir{vector::normalize(ptc.momentum)};
    const point3 glob_pos{sf.local_to_global(ctx, meas.local_position(), dir)};

    // Rough estimate of intersection distance from origin
    const scalar path{vector::norm(glob_pos)};

    // Build an intersection
    using nav_link_t = typename intersection_t::nav_link_type;
    intersection_t intr{
        sf_desc,
        path,
        {meas.local_position()[0], meas.local_position()[1], 0.f},
        static_cast<nav_link_t>(geo_id.volume()),
        detray::intersection::status::e_inside,
        true};

    // TODO: Don't use initial particle momentum
    candidates.emplace_back(glob_pos, dir, intr, q, p);
  }

  // Sort records by intersection distance to origin of the trajectory
  auto sort_by_path = [&](const candidate_type<detector_t>& a,
                          const candidate_type<detector_t>& b) -> bool {
    return (a.intersection < b.intersection);
  };

  std::ranges::stable_sort(candidates, sort_by_path);

  return candidates;
}

}  // namespace traccc::propagation_validator
