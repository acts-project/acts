// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/definitions/units.hpp"

// Detray test include(s)
#include "detray/test/common/build_telescope_detector.hpp"
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
  // Telescope detector configuration
  //
  using metadata_t = test::default_telescope_metadata;
  using tel_detector_t = detector<metadata_t>;
  using test_algebra = metadata_t::algebra_type;
  using scalar = dscalar<test_algebra>;

  /// Set a consistent minimum step size across all tests
  const float min_stepsize{stepping::config{}.min_stepsize};

  tel_det_config<test_algebra, rectangle2D> tel_cfg{20.f * unit<scalar>::mm,
                                                    20.f * unit<scalar>::mm};
  tel_cfg.n_surfaces(10u)
      .length(500.f * unit<scalar>::mm)
      .envelope(500.f * unit<scalar>::um);

  vecmem::host_memory_resource host_mr;

  const auto [tel_det, tel_names] =
      build_telescope_detector<test_algebra>(host_mr, tel_cfg);

  auto white_board = std::make_shared<test::whiteboard>();
  tel_detector_t::geometry_context ctx{};

  // General data consistency of the detector
  test::consistency_check<tel_detector_t>::config cfg_cons{};
  test::register_checks<test::consistency_check>(
      tel_det, tel_names, cfg_cons.name("telescope_detector_consistency"), ctx);

  // Navigation link consistency, discovered by ray intersection
  test::ray_scan<tel_detector_t>::config cfg_ray_scan{};
  cfg_ray_scan.name("telescope_detector_ray_scan");
  cfg_ray_scan.track_generator().n_tracks(10000u);
  cfg_ray_scan.overlaps_tol(min_stepsize);
  // The first surface is at z=0, so shift the track origin back
  cfg_ray_scan.track_generator().origin(0.f, 0.f, -0.05f * unit<scalar>::mm);
  cfg_ray_scan.track_generator().theta_range(0.f,
                                             0.25f * constant<scalar>::pi_4);

  test::register_checks<test::ray_scan>(tel_det, tel_names, cfg_ray_scan, ctx,
                                        white_board);

  // Comparison of straight line navigation with ray scan
  test::straight_line_navigation<tel_detector_t>::config cfg_str_nav{};
  cfg_str_nav.name("telescope_detector_straight_line_navigation");
  cfg_str_nav.n_tracks(cfg_ray_scan.track_generator().n_tracks());
  cfg_str_nav.propagation().stepping.min_stepsize = min_stepsize;
  cfg_str_nav.propagation().navigation.estimate_scattering_noise = false;
  cfg_str_nav.propagation().navigation.intersection.min_mask_tolerance =
      static_cast<float>(cfg_ray_scan.mask_tolerance());
  cfg_str_nav.propagation().navigation.intersection.max_mask_tolerance =
      static_cast<float>(cfg_ray_scan.mask_tolerance());

  test::register_checks<test::straight_line_navigation>(
      tel_det, tel_names, cfg_str_nav, ctx, white_board);

  // Navigation link consistency, discovered by helix intersection
  test::helix_scan<tel_detector_t>::config cfg_hel_scan{};
  cfg_hel_scan.name("telescope_detector_helix_scan");
  // Let the Newton algorithm dynamically choose tol. based on approx. error
  cfg_hel_scan.mask_tolerance(detray::detail::invalid_value<scalar>());
  cfg_hel_scan.track_generator().n_tracks(10000u);
  cfg_hel_scan.overlaps_tol(min_stepsize);
  cfg_hel_scan.track_generator().p_tot(10.f * unit<scalar>::GeV);
  cfg_hel_scan.track_generator().origin(0.f, 0.f, -0.05f * unit<scalar>::mm);
  cfg_hel_scan.track_generator().theta_range(0.f,
                                             0.25f * constant<scalar>::pi_4);

  test::register_checks<test::helix_scan>(tel_det, tel_names, cfg_hel_scan, ctx,
                                          white_board);

  // Comparison of navigation in a constant B-field with helix
  test::helix_navigation<tel_detector_t>::config cfg_hel_nav{};
  cfg_hel_nav.name("telescope_detector_helix_navigation");
  cfg_hel_nav.n_tracks(cfg_hel_scan.track_generator().n_tracks());
  cfg_hel_nav.propagation().stepping.min_stepsize = min_stepsize;
  cfg_hel_nav.propagation().navigation.estimate_scattering_noise = false;
  cfg_hel_nav.propagation().navigation.intersection.overstep_tolerance =
      -100.f * unit<float>::um;

  test::register_checks<test::helix_navigation>(tel_det, tel_names, cfg_hel_nav,
                                                ctx, white_board);

  // Run the material validation
  test::material_scan<tel_detector_t>::config mat_scan_cfg{};
  mat_scan_cfg.name("telescope_detector_material_scan");
  mat_scan_cfg.track_generator().uniform_eta(true).eta_range(1.f, 6.f);
  mat_scan_cfg.track_generator().origin(0.f, 0.f, -0.05f * unit<scalar>::mm);
  mat_scan_cfg.track_generator().phi_steps(100).eta_steps(100);
  mat_scan_cfg.overlaps_tol(min_stepsize);

  // Record the material using a ray scan
  test::register_checks<test::material_scan>(tel_det, tel_names, mat_scan_cfg,
                                             ctx, white_board);

  // Now trace the material during navigation and compare
  test::material_validation<tel_detector_t>::config mat_val_cfg{};
  mat_val_cfg.name("telescope_detector_material_validaiton");
  mat_val_cfg.propagation() = cfg_str_nav.propagation();
  mat_val_cfg.propagation().stepping.min_stepsize = min_stepsize;

  test::register_checks<test::material_validation>(
      tel_det, tel_names, mat_val_cfg, ctx, white_board);

  // Run the checks
  return RUN_ALL_TESTS();
}
