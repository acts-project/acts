// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/cpu/detector_scan.hpp"
#include "detray/test/cpu/material_scan.hpp"
#include "detray/test/device/cuda/material_validation.hpp"
#include "detray/test/device/cuda/navigation_validation.hpp"
#include "detray/test/framework/register_checks.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/framework/whiteboard.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// GTest include
#include <gtest/gtest.h>

// System include(s)
#include <limits>

using namespace detray;

int main(int argc, char **argv) {
  using namespace detray;

  // Filter out the google test flags
  ::testing::InitGoogleTest(&argc, argv);

  using toy_detector_t = detector<test::toy_metadata>;
  using test_algebra = typename toy_detector_t::algebra_type;
  using scalar = dscalar<test_algebra>;

  /// Set a consistent minimum step size across all tests
  const float min_stepsize{stepping::config{}.min_stepsize};

  /// Vecmem memory resource for the device allocations
  vecmem::cuda::device_memory_resource dev_mr{};

  //
  // Toy detector configuration
  //
  toy_det_config<scalar> toy_cfg{};
  toy_cfg.n_brl_layers(4u).n_edc_layers(7u);

  // Build the geometry
  vecmem::host_memory_resource host_mr;
  auto [toy_det, toy_names] =
      build_toy_detector<test_algebra>(host_mr, toy_cfg);

  auto white_board = std::make_shared<test::whiteboard>();
  toy_detector_t::geometry_context ctx{};

  // Navigation link consistency, discovered by ray intersection
  test::ray_scan<toy_detector_t>::config cfg_ray_scan{};
  cfg_ray_scan.name("toy_detector_ray_scan_for_cuda");
  cfg_ray_scan.track_generator().n_tracks(10000u);
  cfg_ray_scan.overlaps_tol(min_stepsize);

  test::register_checks<test::ray_scan>(toy_det, toy_names, cfg_ray_scan, ctx,
                                        white_board);

  // Comparison of straight line navigation with ray scan
  detray::cuda::straight_line_navigation<toy_detector_t>::config cfg_str_nav{};
  cfg_str_nav.name("toy_detector_straight_line_navigation_cuda");
  cfg_str_nav.n_tracks(cfg_ray_scan.track_generator().n_tracks());
  cfg_str_nav.propagation().stepping.min_stepsize = min_stepsize;
  cfg_str_nav.propagation().navigation.estimate_scattering_noise = false;
  cfg_str_nav.propagation().navigation.search_window = {3u, 3u};
  cfg_str_nav.propagation().navigation.intersection.min_mask_tolerance =
      static_cast<float>(cfg_ray_scan.mask_tolerance());
  cfg_str_nav.propagation().navigation.intersection.max_mask_tolerance =
      static_cast<float>(cfg_ray_scan.mask_tolerance());

  test::register_checks<detray::cuda::straight_line_navigation>(
      toy_det, toy_names, cfg_str_nav, ctx, white_board);

  // Navigation link consistency, discovered by helix intersection
  test::helix_scan<toy_detector_t>::config cfg_hel_scan{};
  cfg_hel_scan.name("toy_detector_helix_scan_for_cuda");
  // Let the Newton algorithm dynamically choose tol. based on approx. error
  cfg_hel_scan.mask_tolerance(detray::detail::invalid_value<scalar>());
  // Run only 1000 track in double precision in the CI (time limit)
  if constexpr (std::same_as<scalar, double>) {
    cfg_hel_scan.track_generator().n_tracks(1000u);
  } else {
    cfg_hel_scan.track_generator().n_tracks(10000u);
  }
  cfg_hel_scan.overlaps_tol(min_stepsize);
  cfg_hel_scan.track_generator().randomize_charge(true);
  cfg_hel_scan.track_generator().eta_range(-4.f, 4.f);
  cfg_hel_scan.track_generator().p_T(1.f * unit<scalar>::GeV);

  test::register_checks<test::helix_scan>(toy_det, toy_names, cfg_hel_scan, ctx,
                                          white_board);

  // Comparison of navigation in a constant B-field with helix
  detray::cuda::helix_navigation<toy_detector_t>::config cfg_hel_nav{};
  cfg_hel_nav.name("toy_detector_helix_navigation_cuda");
  cfg_hel_nav.n_tracks(cfg_hel_scan.track_generator().n_tracks());
  cfg_hel_nav.propagation().stepping.min_stepsize = min_stepsize;
  cfg_hel_nav.propagation().navigation.estimate_scattering_noise = false;
  cfg_hel_nav.propagation().navigation.search_window = {3u, 3u};

  test::register_checks<detray::cuda::helix_navigation>(
      toy_det, toy_names, cfg_hel_nav, ctx, white_board);

  // Run the material validation - Material Maps
  test::material_scan<toy_detector_t>::config mat_scan_cfg{};
  mat_scan_cfg.name("toy_detector_material_scan_for_cuda");
  mat_scan_cfg.track_generator().uniform_eta(true).eta_range(-4.f, 4.f);
  mat_scan_cfg.track_generator().phi_steps(100).eta_steps(100);
  mat_scan_cfg.overlaps_tol(min_stepsize);

  // Record the material using a ray scan
  test::register_checks<test::material_scan>(toy_det, toy_names, mat_scan_cfg,
                                             ctx, white_board);

  // Now trace the material during navigation and compare
  detray::cuda::material_validation<toy_detector_t>::config mat_val_cfg{};
  mat_val_cfg.name("toy_detector_material_validaiton_cuda");
  mat_val_cfg.device_mr(&dev_mr);
  mat_val_cfg.tol(1e-5f);  // < Reduce tolerance for single precision tests
  mat_val_cfg.propagation() = cfg_str_nav.propagation();
  mat_val_cfg.propagation().stepping.min_stepsize = min_stepsize;

  test::register_checks<detray::cuda::material_validation>(
      toy_det, toy_names, mat_val_cfg, ctx, white_board);

  // Run the material validation - Homogeneous material
  toy_cfg.use_material_maps(false);

  auto [toy_det_hom_mat, toy_names_hom_mat] =
      build_toy_detector<test_algebra>(host_mr, toy_cfg);
  toy_names_hom_mat.set_detector_name(toy_names_hom_mat.get_detector_name() +
                                      "_hom_material");

  // Record the material using a ray scan
  mat_scan_cfg.name("toy_detector_hom_material_scan_for_cuda");
  test::register_checks<test::material_scan>(toy_det_hom_mat, toy_names_hom_mat,
                                             mat_scan_cfg, ctx, white_board);

  // Now trace the material during navigation and compare
  mat_val_cfg.name("toy_detector_hom_material_validaiton_cuda");
  test::register_checks<detray::cuda::material_validation>(
      toy_det_hom_mat, toy_names_hom_mat, mat_val_cfg, ctx, white_board);

  // Run the checks
  return RUN_ALL_TESTS();
}
