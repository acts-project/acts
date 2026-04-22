// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"

// Detray test include(s)
#include "detray/test/common/build_wire_chamber.hpp"
#include "detray/test/cpu/detector_consistency.hpp"
#include "detray/test/cpu/detector_scan.hpp"
#include "detray/test/cpu/material_scan.hpp"
#include "detray/test/cpu/material_validation.hpp"
#include "detray/test/cpu/navigation_validation.hpp"
#include "detray/test/framework/register_checks.hpp"
#include "detray/test/framework/whiteboard.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

int main(int argc, char **argv) {
  using namespace detray;

  // Filter out the google test flags
  ::testing::InitGoogleTest(&argc, argv);

  //
  // Wire Chamber configuration
  //
  vecmem::host_memory_resource host_mr;

  using metadata_t = test::wire_chamber_metadata;
  using wire_chamber_t = detector<metadata_t>;
  using test_algebra = metadata_t::algebra_type;
  using scalar = dscalar<test_algebra>;

  /// Set a consistent minimum step size across all tests
  const float min_stepsize{stepping::config{}.min_stepsize};

  wire_chamber_config<scalar> wire_chamber_cfg{};
  wire_chamber_cfg.half_z(500.f * unit<scalar>::mm);

  std::clog << wire_chamber_cfg << std::endl;

  auto [det, names] =
      build_wire_chamber<test_algebra>(host_mr, wire_chamber_cfg);

  auto white_board = std::make_shared<test::whiteboard>();
  wire_chamber_t::geometry_context ctx{};

  // General data consistency of the detector
  test::consistency_check<wire_chamber_t>::config cfg_cons{};
  test::register_checks<test::consistency_check>(
      det, names, cfg_cons.name("wire_chamber_consistency"), ctx);

  // Navigation link consistency, discovered by ray intersection
  test::ray_scan<wire_chamber_t>::config cfg_ray_scan{};
  cfg_ray_scan.name("wire_chamber_ray_scan");
  cfg_ray_scan.track_generator().seed(42u);
  cfg_ray_scan.track_generator().n_tracks(10000u);
  cfg_ray_scan.overlaps_tol(min_stepsize);

  test::register_checks<test::ray_scan>(det, names, cfg_ray_scan, ctx,
                                        white_board);

  // Comparison of straight line navigation with ray scan
  test::straight_line_navigation<wire_chamber_t>::config cfg_str_nav{};
  cfg_str_nav.name("wire_chamber_straight_line_navigation");
  cfg_str_nav.n_tracks(cfg_ray_scan.track_generator().n_tracks());
  cfg_str_nav.propagation().stepping.min_stepsize = min_stepsize;
  cfg_str_nav.propagation().navigation.estimate_scattering_noise = false;
  cfg_str_nav.propagation().navigation.search_window = {3u, 3u};
  cfg_str_nav.propagation().navigation.intersection.min_mask_tolerance =
      static_cast<float>(cfg_ray_scan.mask_tolerance());
  cfg_str_nav.propagation().navigation.intersection.max_mask_tolerance =
      static_cast<float>(cfg_ray_scan.mask_tolerance());

  test::register_checks<test::straight_line_navigation>(det, names, cfg_str_nav,
                                                        ctx, white_board);

  // Navigation link consistency, discovered by helix intersection
  test::helix_scan<wire_chamber_t>::config cfg_hel_scan{};
  cfg_hel_scan.name("wire_chamber_helix_scan");
  // Let the Newton algorithm dynamically choose tol. based on approx. error
  cfg_hel_scan.mask_tolerance(detray::detail::invalid_value<scalar>());
  cfg_hel_scan.track_generator().n_tracks(10000u);
  cfg_hel_scan.track_generator().randomize_charge(true);
  cfg_hel_scan.track_generator().eta_range(-1.f, 1.f);
  // TODO: Fails for smaller momenta
  cfg_hel_scan.track_generator().p_T(4.f * unit<scalar>::GeV);
  cfg_hel_scan.overlaps_tol(min_stepsize);

  test::register_checks<test::helix_scan>(det, names, cfg_hel_scan, ctx,
                                          white_board);

  // Comparison of navigation in a constant B-field with helix
  test::helix_navigation<wire_chamber_t>::config cfg_hel_nav{};
  cfg_hel_nav.name("wire_chamber_helix_navigation");
  cfg_hel_nav.n_tracks(cfg_hel_scan.track_generator().n_tracks());
  cfg_hel_nav.propagation().stepping.min_stepsize = min_stepsize;
  cfg_hel_nav.propagation().navigation.estimate_scattering_noise = false;
  cfg_hel_nav.propagation().navigation.intersection.min_mask_tolerance *= 11.f;
  cfg_hel_nav.propagation().navigation.search_window = {3u, 3u};

  test::register_checks<test::helix_navigation>(det, names, cfg_hel_nav, ctx,
                                                white_board);

  // Run the material validation
  test::material_scan<wire_chamber_t>::config mat_scan_cfg{};
  mat_scan_cfg.name("wire_chamber_material_scan");
  mat_scan_cfg.track_generator().uniform_eta(true).eta_range(-1.f, 1.f);
  mat_scan_cfg.track_generator().phi_steps(100).eta_steps(100);
  mat_scan_cfg.overlaps_tol(min_stepsize);

  // Record the material using a ray scan
  test::register_checks<test::material_scan>(det, names, mat_scan_cfg, ctx,
                                             white_board);

  // Now trace the material during navigation and compare
  test::material_validation<wire_chamber_t>::config mat_val_cfg{};
  mat_val_cfg.name("wire_chamber_material_validaiton");
  // Reduce tolerance for single precision tests
  if constexpr (std::is_same_v<scalar, float>) {
    mat_val_cfg.relative_error(130.f);
  }
  mat_val_cfg.propagation() = cfg_str_nav.propagation();
  mat_val_cfg.propagation().stepping.min_stepsize = min_stepsize;

  test::register_checks<test::material_validation>(det, names, mat_val_cfg, ctx,
                                                   white_board);

  // Run the checks
  return RUN_ALL_TESTS();
}
