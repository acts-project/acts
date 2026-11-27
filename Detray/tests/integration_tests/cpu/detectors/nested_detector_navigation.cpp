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
#include "detray/test/common/build_nested_detector.hpp"
#include "detray/test/cpu/detector_consistency.hpp"
#include "detray/test/cpu/detector_scan.hpp"
#include "detray/test/cpu/navigation_validation.hpp"
#include "detray/test/framework/register_checks.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/framework/whiteboard.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <string>

int main(int argc, char **argv) {
  using namespace detray;

  // Filter out the google test flags
  ::testing::InitGoogleTest(&argc, argv);

  //
  // Nested detector configuration
  //
  using test_algebra = test::algebra;
  using scalar = test::scalar;

  /// Set a consistent minimum step size across all tests
  const float min_stepsize{stepping::config{}.min_stepsize};

  vecmem::host_memory_resource host_mr;

  // Daughter volumes that are flush with their mother boundary
  const auto [nested_det, nested_names] =
      build_nested_detector<test_algebra>(host_mr);

  using detector_t = std::remove_cvref_t<decltype(nested_det)>;

  auto white_board = std::make_shared<test::whiteboard>();
  typename detector_t::geometry_context ctx{};

  // The track origin is in the world volume (but not on a surface)
  const scalar ori_x{0.37f * unit<scalar>::mm};
  const scalar ori_y{-0.61f * unit<scalar>::mm};
  const scalar ori_z{-20.3f * unit<scalar>::mm};

  // General data consistency of the detector
  test::consistency_check<detector_t>::config cfg_cons{};
  test::register_checks<test::consistency_check>(
      nested_det, nested_names, cfg_cons.name("nested_detector_consistency"),
      ctx);

  // Navigation link consistency, discovered by ray intersection
  test::ray_scan<detector_t>::config cfg_ray_scan{};
  cfg_ray_scan.name("nested_detector_ray_scan");
  cfg_ray_scan.track_generator().n_tracks(10000u);
  cfg_ray_scan.track_generator().origin(ori_x, ori_y, ori_z);
  cfg_ray_scan.overlaps_tol(min_stepsize);
  // Portals of mother and daughter volumes can be at the same position
  cfg_ray_scan.overlaps_removal(false);

  test::register_checks<test::ray_scan>(nested_det, nested_names, cfg_ray_scan,
                                        ctx, white_board);

  // Comparison of straight line navigation with ray scan: The ray scan
  // also finds the portals of mother volumes, which are covered by daughter
  // volumes. Compare only the sensitive surfaces.
  test::straight_line_navigation<detector_t>::config cfg_str_nav{};
  cfg_str_nav.name("nested_detector_straight_line_navigation");
  cfg_str_nav.n_tracks(cfg_ray_scan.track_generator().n_tracks());
  cfg_str_nav.collect_sensitives_only(true);
  cfg_str_nav.propagation().stepping.min_stepsize = min_stepsize;
  cfg_str_nav.propagation().navigation.estimate_scattering_noise = false;
  cfg_str_nav.propagation().navigation.intersection.min_mask_tolerance =
      static_cast<float>(cfg_ray_scan.mask_tolerance());
  cfg_str_nav.propagation().navigation.intersection.max_mask_tolerance =
      static_cast<float>(cfg_ray_scan.mask_tolerance());

  test::register_checks<test::straight_line_navigation>(
      nested_det, nested_names, cfg_str_nav, ctx, white_board);

  // Navigation link consistency, discovered by helix intersection
  test::helix_scan<detector_t>::config cfg_hel_scan{};
  cfg_hel_scan.name("nested_detector_helix_scan");
  // Let the Newton algorithm dynamically choose tol. based on approx. error
  cfg_hel_scan.mask_tolerance(detray::detail::invalid_value<scalar>());
  cfg_hel_scan.track_generator().n_tracks(10000u);
  cfg_hel_scan.track_generator().p_tot(10.f * unit<scalar>::GeV);
  cfg_hel_scan.track_generator().origin(ori_x, ori_y, ori_z);
  cfg_hel_scan.overlaps_tol(min_stepsize);
  cfg_hel_scan.overlaps_removal(false);

  test::register_checks<test::helix_scan>(nested_det, nested_names,
                                          cfg_hel_scan, ctx, white_board);

  // Comparison of navigation in a constant B-field with helix
  test::helix_navigation<detector_t>::config cfg_hel_nav{};
  cfg_hel_nav.name("nested_detector_helix_navigation");
  cfg_hel_nav.n_tracks(cfg_hel_scan.track_generator().n_tracks());
  cfg_hel_nav.collect_sensitives_only(true);
  cfg_hel_nav.propagation().stepping.min_stepsize = min_stepsize;
  cfg_hel_nav.propagation().navigation.estimate_scattering_noise = false;
  cfg_hel_nav.propagation().navigation.intersection.overstep_tolerance =
      -100.f * unit<float>::um;

  test::register_checks<test::helix_navigation>(nested_det, nested_names,
                                                cfg_hel_nav, ctx, white_board);

  // Run the checks
  return RUN_ALL_TESTS();
}
