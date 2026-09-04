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
#include "detray/utils/logging.hpp"

// Detray plugin include(s)
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// Detray IO include(s)
#include "detray/io/utils/create_path.hpp"

// System include(s)
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <unordered_set>

namespace detray::detail {

/// Find the unique volume indices that the trajectory crossed
/// - intersection record
template <typename record_t>
std::unordered_set<dindex> get_volume_indices(
    const dvector<record_t> &intersection_record) {
  std::unordered_set<dindex> volumes{};
  volumes.reserve(intersection_record.size());
  for (const auto &single_ir : intersection_record) {
    if (single_ir.intersection.surface().identifier().is_invalid()) {
      continue;
    }
    volumes.insert(single_ir.vol_idx);
  }

  return volumes;
}

/// Find the unique volume indices that the trajectory crossed
/// - intersection collection
template <typename surface_t, concepts::algebra algebra_t>
std::unordered_set<dindex> get_volume_indices(
    const dvector<detray::intersection2D<
        surface_t, algebra_t, intersection::contains_pos>> &intersections) {
  std::unordered_set<dindex> volumes{};
  volumes.reserve(intersections.size());
  for (const auto &intr : intersections) {
    if (intr.surface().identifier().is_invalid()) {
      continue;
    }
    volumes.insert(intr.surface().volume());
  }

  return volumes;
}

/// Find the unique volume indices that the trajectory crossed
/// - object tracer candidate collection
template <template <typename> class candidate_t, typename surface_t,
          concepts::algebra algebra_t>
std::unordered_set<dindex> get_volume_indices(
    const dvector<candidate_t<detray::intersection2D<
        surface_t, algebra_t, intersection::contains_pos>>> &candidates) {
  std::unordered_set<dindex> volumes{};
  volumes.reserve(candidates.size());
  for (const auto &cand : candidates) {
    if (cand.intersection.surface().identifier().is_invalid()) {
      continue;
    }
    volumes.insert(cand.intersection.surface().volume());
  }

  return volumes;
}

/// Transcribe the intersections from an intersection trace into a standalone
/// vector.
///
/// @param intersection_trace the input intersection trace
///
/// @returns a vector of intersections in the same order as the input trace.
template <typename record_t, typename ALLOC>
auto transcribe_intersections(
    const std::vector<record_t, ALLOC> &intersection_trace) {
  using intersection_t = typename record_t::intersection_type;

  std::vector<intersection_t> intersections{};
  intersections.reserve(intersection_trace.size());
  for (auto &ir : intersection_trace) {
    // Dummy record
    if (ir.intersection.surface().identifier().is_invalid()) {
      continue;
    }
    intersections.push_back(ir.intersection);
  }

  return intersections;
}

/// @returns the svg of the intersections (truth and track) and the trajectory
template <typename detector_t, typename truth_trace_t, class traj_t,
          typename recorded_trace_t, typename view_t>
auto draw_intersection_and_traj_svg(
    const typename detector_t::geometry_context gctx,
    detray::svgtools::illustrator<detector_t> &il,
    const truth_trace_t &truth_trace, const traj_t &traj,
    const std::string &traj_name, const recorded_trace_t &recorded_trace,
    const view_t &view,
    dindex_range highlight_idx = {detail::invalid_value<dindex>(),
                                  detail::invalid_value<dindex>()}) {
  // Get only the intersections from the traces
  auto truth_intersections = transcribe_intersections(truth_trace);
  auto recorded_intersections = transcribe_intersections(recorded_trace);

  // Was a dummy record removed before printing?
  if (truth_intersections.size() < truth_trace.size()) {
    highlight_idx.shift(-1);
  }

  // Draw the truth intersections as black crosses
  actsvg::svg::object svg_traj;
  if (!truth_intersections.empty()) {
    svg_traj = il.draw_intersections("truth_trace", truth_intersections,
                                     traj.dir(0.f), view, highlight_idx, gctx);
  }

  // Draw an approximation of the trajectory with the recorded intersections
  const auto path{math::max(math::fabs(truth_intersections.front().path()),
                            math::fabs(truth_intersections.back().path()))};
  if (!recorded_intersections.empty()) {
    svg_traj.add_object(il.draw_intersections_and_trajectory(
        traj_name, recorded_intersections, traj, view, path, highlight_idx,
        gctx));
  } else {
    svg_traj.add_object(il.draw_trajectory(traj_name, traj, path, view));
  }

  return svg_traj;
}

/// Display the geometry, intersection and track data via @c svgtools
template <typename detector_t, typename truth_trace_t, class traj_t,
          typename recorded_trace_t>
inline void svg_display(
    const typename detector_t::geometry_context gctx,
    detray::svgtools::illustrator<detector_t> &il,
    const truth_trace_t &truth_trace, const traj_t &traj,
    const std::string &traj_name,
    const std::string &outfile = "detector_display",
    const recorded_trace_t &recorded_trace = {},
    dindex_range highlight_idx = {detail::invalid_value<dindex>(),
                                  detail::invalid_value<dindex>()},
    const bool verbose = true, const std::string &outdir = "./plots/") {
  // Gather all volumes that need to be displayed
  auto volumes = get_volume_indices(truth_trace);
  if (!recorded_trace.empty()) {
    const auto more_volumes = get_volume_indices(truth_trace);
    volumes.insert(more_volumes.begin(), more_volumes.end());
  }

  // General options
  auto path = detray::io::create_path(outdir);

  actsvg::style::stroke stroke_black = actsvg::style::stroke();

  // x-y axis.
  auto xy_axis = actsvg::draw::x_y_axes("axes", {-1100, 1100}, {-1100, 1100},
                                        stroke_black, "x", "y");
  // z-r axis.
  auto zr_axis = actsvg::draw::x_y_axes("axes", {-3100, 3100}, {-5, 1100},
                                        stroke_black, "z", "r");
  // Creating the views.
  const actsvg::views::x_y xy{};
  const actsvg::views::z_r zr{};

  // xy - view
  auto svg_traj =
      draw_intersection_and_traj_svg(gctx, il, truth_trace, traj, traj_name,
                                     recorded_trace, xy, highlight_idx);

  const auto [vol_xy_svg, _] = il.draw_volumes(volumes, xy, gctx);
  detray::svgtools::write_svg(
      path / (outfile + "_" + vol_xy_svg._id + "_" + traj_name),
      {xy_axis, vol_xy_svg, svg_traj});

  // zr - view
  svg_traj =
      draw_intersection_and_traj_svg(gctx, il, truth_trace, traj, traj_name,
                                     recorded_trace, zr, highlight_idx);

  const auto vol_zr_svg = il.draw_detector(zr, gctx);
  detray::svgtools::write_svg(
      path / (outfile + "_" + vol_zr_svg._id + "_" + traj_name),
      {zr_axis, vol_zr_svg, svg_traj});

  if (verbose) {
    DETRAY_INFO_HOST("Wrote svgs for debugging in: " << path << "\n");
  }
}

}  // namespace detray::detail
