// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/navigation/volume_graph.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/logging.hpp"

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/track_generators.hpp"
#include "detray/test/utils/hash_tree.hpp"
#include "detray/test/validation/detector_scan_utils.hpp"
#include "detray/test/validation/detector_scanner.hpp"

// Example linear algebra plugin: std::array
#include "detray/tutorial/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>

// Hash of the "correct" geometry
constexpr std::size_t root_hash = 11359580520962287982ul;

/// Check a given detecor for consistent linking by shooting rays/helices and
/// recording every intersection with the geometry. This intersection record
/// can then be checked for matching portals at the volume boundary surfaces
/// ( @c trace_intersections ) and checked for a consistent 'path' from volume
/// to volume ( @c check_consistency ). See also documentation in
/// 'tests/common/tools/ray_detector_scan_utils.hpp'.
int main() {
  using algebra_t = detray::tutorial::algebra_t;
  using scalar = detray::tutorial::scalar;
  using detector_t = detray::detector<detray::toy_metadata<algebra_t>>;

  // Can also be performed with helices
  using ray_t = detray::detail::ray<algebra_t>;

  std::clog << "Ray Scan Tutorial\n=================\n\n";

  // Build the geometry
  vecmem::host_memory_resource host_mr;
  const auto [det, names] = detray::build_toy_detector<algebra_t>(host_mr);

  // The current geometry context
  const detector_t::geometry_context gctx{};

  // Visualization style to be applied to the SVGs
  detray::svgtools::styling::style svg_style =
      detray::svgtools::styling::tableau_colorblind::style;

  // Optional: get the volume adjaceny matrix from ray scan
  detray::volume_graph graph(det);
  const auto &adj_mat = graph.adjacency_matrix();  // < need this for the size
  detray::dvector<detray::dindex> adj_mat_scan(adj_mat.size(), 0);

  // Keep track of the objects that have already been seen per volume
  std::unordered_set<detray::dindex> obj_hashes = {};

  // Index of the volume that the ray origin lies in
  detray::dindex start_index{0u};
  std::size_t n_rays{0u};

  // Generate a number of random rays
  using generator_t =
      detray::detail::random_numbers<scalar,
                                     std::uniform_real_distribution<scalar>>;
  using ray_generator_t = detray::random_track_generator<ray_t, generator_t>;

  ray_generator_t::configuration cfg{};
  cfg.n_tracks(10000).p_T(1.f * detray::unit<scalar>::GeV);

  ray_generator_t ray_generator{cfg};

  // Run the check
  DETRAY_INFO_HOST("\nScanning " << det.name(names) << " ("
                                 << ray_generator.size() << " rays) ...\n");
  using intersection_trace_t =
      detray::dvector<detray::intersection_record<detector_t>>;

  bool success = true;
  for (const auto ray : ray_generator) {
    // Record all intersections and surfaces along the ray
    const auto intersection_trace =
        detray::detector_scanner::run<detray::ray_scan>(gctx, det, ray);

    bool check_result = detray::detector_scanner::check_trace<detector_t>(
        intersection_trace, start_index, adj_mat_scan, obj_hashes);

    if (!check_result) {
      detray::detector_scanner::display_error(
          gctx, det, names, "ray_scan_tutorial", ray, intersection_trace,
          svg_style, n_rays, ray_generator.size(), intersection_trace_t{});
    }
    success = success && check_result;

    ++n_rays;
  }

  // Check result
  DETRAY_INFO_HOST("Ray scan: " << (success ? "OK" : "FAILURE"));

  // Compare the adjacency that was discovered in the ray scan to the hashed
  // one for the toy detector.
  // The hash tree is still Work in Progress !
  auto geo_checker = detray::hash_tree(adj_mat);
  const bool check_links = (geo_checker.root() == root_hash);

  DETRAY_INFO_HOST("All links reachable: " << (check_links ? "OK" : "FAILURE"));
}
