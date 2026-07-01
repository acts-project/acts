// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/navigation/detail/intersection_kernel.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersector.hpp"
#include "detray/tracks/free_track_parameters.hpp"
#include "detray/tracks/trajectories.hpp"

// Detray IO include(s)
#include "detray/io/csv/intersection2D.hpp"
#include "detray/io/csv/track_parameters.hpp"

// Test include(s)
#include "detray/test/utils/data_record.hpp"

// System include(s)
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace detray {

/// Record of a surface intersection along a track
/*template <typename detector_t>
struct intersection_record {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using track_parameter_type = free_track_parameters<algebra_t>;
  using intersection_type =
      intersection2D<typename detector_t::surface_type, algebra_t,
                     intersection::contains_pos>;

  /// The charge associated with the track parameters
  scalar_t charge{};
  /// Current global track parameters
  track_parameter_type track_param{
      {0.f, 0.f, 0.f}, 0.f, {0.f, 0.f, 1.f}, detail::invalid_value<scalar_t>()};
  /// Index of the volume the intersection was found in
  dindex vol_idx{};
  /// The intersection result, including the surface descriptor
  intersection_type intersection{};
};*/

/// @brief struct that holds functionality to shoot a parametrized particle
/// trajectory through a detector.
///
/// Records intersections with every detector surface along the trajectory.
template <typename trajectory_t>
struct brute_force_scan {
  template <typename D>
  using intersection_trace_type = std::vector<intersection_record<D>>;
  using trajectory_type = trajectory_t;

  template <typename detector_t>
  inline auto operator()(const typename detector_t::geometry_context ctx,
                         const detector_t &detector, const trajectory_t &traj,
                         const typename detector_t::scalar_type mask_tol = 0.f,
                         const typename detector_t::scalar_type p =
                             1.f *
                             unit<typename detector_t::scalar_type>::GeV) {
    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    using sf_desc_t = typename detector_t::surface_type;
    using nav_link_t = typename detector_t::surface_type::navigation_link;

    using intersection_t =
        typename intersection_record<detector_t>::intersection_type;
    using intersection_kernel_t = detail::intersection_initialize<intersector>;

    constexpr scalar_t external_mask_tol{0.f};
    const intersection::config intr_cfg{
        .min_mask_tolerance = static_cast<float>(mask_tol),
        .max_mask_tolerance = static_cast<float>(mask_tol),
        .mask_tolerance_scalor = 0.f,
        .overstep_tolerance = 0.f};

    intersection_trace_type<detector_t> intersection_trace;

    const auto &trf_store = detector.transform_store();

    assert(p > 0.f);
    const scalar_t q{p * traj.qop()};

    std::vector<intersection_t> intersections{};
    intersections.reserve(100u);

    // Loop over all surfaces in the detector
    for (const sf_desc_t &sf_desc : detector.surfaces()) {
      // Retrieve candidate(s) from the surface
      const auto sf = geometry::surface{detector, sf_desc};
      sf.template visit_mask<intersection_kernel_t>(
          intersections, traj, sf_desc, trf_store, ctx, intr_cfg,
          external_mask_tol);

      // Candidate is invalid if it lies in the opposite direction
      for (auto &sfi : intersections) {
        if (sfi.is_along()) {
          sfi.surface() = sf_desc;
          // Record the intersection
          intersection_trace.emplace_back(traj.pos(sfi.path()),
                                          traj.dir(sfi.path()), sfi, q, p,
                                          sf.volume());
        }
      }
      intersections.clear();
    }

    // Need to have at least an exit portal
    if (intersection_trace.empty()) {
      std::stringstream stream;
      stream << "No intersections found for traj: " << traj << std::endl;
      throw std::runtime_error(stream.str());
    }

    // Save initial track position as dummy intersection record
    const auto &first_record = intersection_trace.front();
    intersection_t start_intersection{};
    // Must not be invalid, since it will otherwise throw the navigation
    // validation off
    start_intersection.set_surface(first_record.intersection.surface());
    start_intersection.surface().set_id(surface_id::e_passive);
    start_intersection.surface().set_index(0);
    start_intersection.surface()
        .material()
        .set_id(detector_t::material::id::e_none)
        .set_index(dindex_invalid);
    start_intersection.set_path(0.f);
    start_intersection.set_local({0.f, 0.f, 0.f});
    start_intersection.set_volume_link(
        static_cast<nav_link_t>(first_record.vol_idx));

    intersection_trace.insert(intersection_trace.begin(),
                              intersection_record<detector_t>{
                                  traj.pos(), traj.dir(), start_intersection, q,
                                  p, first_record.vol_idx});

    return intersection_trace;
  }
};

template <concepts::algebra algebra_t>
using ray_scan = brute_force_scan<detail::ray<algebra_t>>;

template <concepts::algebra algebra_t>
using helix_scan = brute_force_scan<detail::helix<algebra_t>>;

/// Run a scan on detector object by shooting test particles through it
namespace detector_scanner {

template <template <typename> class scan_type, typename detector_t,
          typename trajectory_t, typename... Args>
inline auto run(const typename detector_t::geometry_context gctx,
                const detector_t &detector, const trajectory_t &traj,
                Args &&...args) {
  using algebra_t = typename detector_t::algebra_type;
  using nav_link_t = typename detector_t::surface_type::navigation_link;

  auto intersection_record =
      scan_type<algebra_t>{}(gctx, detector, traj, std::forward<Args>(args)...);

  using record_t = typename decltype(intersection_record)::value_type;

  // HACK: For whatever reason, std::stable_sort really dislikes custom
  // aligned types like the ones in Eigen and Fastor, so we have to sort
  // by indices and then reconstruct the sorted intersection record.
  auto sort_path = [&](const record_t &a, const record_t &b) -> bool {
    return (a.intersection < b.intersection);
  };

  std::ranges::stable_sort(intersection_record, sort_path);

  // Make sure the intersection record terminates at world portals
  auto is_world_exit = [](const record_t &r) {
    return r.intersection.volume_link() ==
           detray::detail::invalid_value<nav_link_t>();
  };

  if (auto it = std::ranges::find_if(intersection_record, is_world_exit);
      it != intersection_record.end()) {
    auto n{static_cast<std::size_t>(it - intersection_record.begin())};
    intersection_record.resize(n + 1u);
  }

  return intersection_record;
}

/// Write the @param intersection_traces to file
template <typename detector_t>
inline auto write_intersections(
    const std::string &intersection_file_name,
    const std::vector<std::vector<intersection_record<detector_t>>>
        &intersection_traces) {
  using record_t = intersection_record<detector_t>;
  using intersection_t = typename record_t::intersection_type;

  std::vector<std::vector<intersection_t>> intersections{};

  // Split data
  for (const auto &trace : intersection_traces) {
    auto &intrs = intersections.emplace_back();
    intrs.reserve(trace.size());

    for (const auto &record : trace) {
      intrs.push_back(record.intersection);
    }
  }

  // Write to file
  io::csv::write_intersection2D(intersection_file_name, intersections);
}

/// Write the @param intersection_traces to file
template <typename record_t>
inline auto write_intersections(
    const std::string &intersection_file_name,
    const dvector<dvector<record_t>> &intersection_traces) {
  using intersection_t = typename record_t::intersection_type;

  std::vector<std::vector<intersection_t>> intersections{};

  // Split data
  for (const auto &trace : intersection_traces) {
    auto &intrs = intersections.emplace_back();
    intrs.reserve(trace.size());

    for (const auto &record : trace) {
      intrs.push_back(record.intersection);
    }
  }

  // Write to file
  io::csv::write_intersection2D(intersection_file_name, intersections);
}

/// Write the @param intersection_traces to file
template <typename detector_t>
inline auto write_tracks(
    const std::string &track_param_file_name,
    const std::vector<std::vector<intersection_record<detector_t>>>
        &intersection_traces) {
  using scalar_t = dscalar<typename detector_t::algebra_type>;
  using track_param_t =
      free_track_parameters<typename detector_t::algebra_type>;

  std::vector<std::vector<std::pair<scalar_t, track_param_t>>> track_params{};

  // Split data
  for (const auto &trace : intersection_traces) {
    track_params.push_back({});
    track_params.back().reserve(trace.size());

    for (const auto &record : trace) {
      track_params.back().emplace_back(record.charge, record.track_param());
    }
  }

  // Write to file
  io::csv::write_free_track_params(track_param_file_name, track_params);
}

/// Read the @param intersection_record from file
template <typename detector_t>
inline auto read(const std::string &intersection_file_name,
                 const std::string &track_param_file_name,
                 std::vector<std::vector<intersection_record<detector_t>>>
                     &intersection_traces) {
  using track_t = free_track_parameters<typename detector_t::algebra_type>;

  // Read from file
  auto intersections_per_track =
      io::csv::read_intersection2D<detector_t>(intersection_file_name);
  auto track_params_per_track =
      io::csv::read_free_track_params<detector_t>(track_param_file_name);

  if (intersections_per_track.size() != track_params_per_track.size()) {
    throw std::invalid_argument(
        "Detector scanner: intersection and track parameters "
        "collections "
        "have different size");
  }

  // Interleave data
  for (dindex trk_idx = 0u; trk_idx < intersections_per_track.size();
       ++trk_idx) {
    const auto &intersections = intersections_per_track[trk_idx];
    const auto &track_params = track_params_per_track[trk_idx];

    // Check track id
    if (intersections.size() != track_params.size()) {
      throw std::invalid_argument(
          "Detector scanner: Found different number of intersections "
          "and "
          "track parameters for track no." +
          std::to_string(trk_idx));
    }

    // Check for empty input traces
    if (intersections.empty()) {
      throw std::invalid_argument("Detector scanner: Found empty trace no." +
                                  std::to_string(trk_idx));
    }

    // Add new trace
    if (intersection_traces.size() <= trk_idx) {
      intersection_traces.push_back({});
    }

    // Add records to trace
    for (dindex i = 0u; i < intersections.size(); ++i) {
      const track_t &free_param = track_params[i].second;
      intersection_traces[trk_idx].push_back(intersection_record<detector_t>{
          free_param.pos(), free_param.dir(), intersections[i],
          track_params[i].first, free_param.p(track_params[i].first),
          intersections[i].surface().volume()});
    }
  }
}

}  // namespace detector_scanner

}  // namespace detray
