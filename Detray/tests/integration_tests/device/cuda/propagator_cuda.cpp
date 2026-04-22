// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "propagator_cuda_kernel.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// GTest include
#include <gtest/gtest.h>

using namespace detray;

class CudaPropConstBFieldMng
    : public ::testing::TestWithParam<std::tuple<float, float, vector3>> {};

/// Propagation test using unified memory
TEST_P(CudaPropConstBFieldMng, propagator) {
  // VecMem memory resource(s)
  vecmem::cuda::managed_memory_resource mng_mr;

  // Test configuration
  propagator_test_config cfg{};
  cfg.track_generator.phi_steps(20).theta_steps(20);
  cfg.track_generator.p_tot(10.f * unit<scalar>::GeV);
  cfg.track_generator.eta_range(-3.f, 3.f);
  cfg.propagation.navigation.search_window = {3u, 3u};
  // Configuration for non-z-aligned B-fields
  cfg.propagation.navigation.intersection.overstep_tolerance =
      std::get<0>(GetParam());
  cfg.propagation.stepping.step_constraint = std::get<1>(GetParam());

  // Get the magnetic field
  const vector3 B = std::get<2>(GetParam());
  auto field = create_const_field<scalar>(B);

  // Create the toy geometry
  auto [det, names] = build_toy_detector<test_algebra>(mng_mr);

  run_propagation_test<bfield::const_bknd_t<scalar>>(
      &mng_mr, det, cfg, detray::get_data(det), std::move(field));
}

class CudaPropConstBFieldCpy
    : public ::testing::TestWithParam<std::tuple<float, float, vector3>> {};

/// Propagation test using vecmem copy
TEST_P(CudaPropConstBFieldCpy, propagator) {
  // VecMem memory resource(s)
  vecmem::host_memory_resource host_mr;
  vecmem::cuda::managed_memory_resource mng_mr;
  vecmem::cuda::device_memory_resource dev_mr;

  vecmem::cuda::copy cuda_cpy;

  // Test configuration
  propagator_test_config cfg{};
  cfg.track_generator.phi_steps(20u).theta_steps(20u);
  cfg.track_generator.p_tot(10.f * unit<scalar>::GeV);
  cfg.track_generator.eta_range(-3.f, 3.f);
  cfg.propagation.navigation.search_window = {3u, 3u};
  // Configuration for non-z-aligned B-fields
  cfg.propagation.navigation.intersection.overstep_tolerance =
      std::get<0>(GetParam());
  cfg.propagation.stepping.step_constraint = std::get<1>(GetParam());

  // Get the magnetic field
  const vector3 B = std::get<2>(GetParam());
  auto field = create_const_field<scalar>(B);

  // Create the toy geometry
  auto [det, names] = build_toy_detector<test_algebra>(host_mr);

  auto det_buff = detray::get_buffer(det, dev_mr, cuda_cpy);

  run_propagation_test<bfield::const_bknd_t<scalar>>(
      &mng_mr, det, cfg, detray::get_data(det_buff), std::move(field));
}

INSTANTIATE_TEST_SUITE_P(
    CudaPropagatorValidation1, CudaPropConstBFieldMng,
    ::testing::Values(std::make_tuple(-100.f * unit<float>::um,
                                      std::numeric_limits<float>::max(),
                                      vector3{0.f * unit<scalar>::T,
                                              0.f * unit<scalar>::T,
                                              2.f * unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    CudaPropagatorValidation2, CudaPropConstBFieldMng,
    ::testing::Values(std::make_tuple(-400.f * unit<float>::um,
                                      std::numeric_limits<float>::max(),
                                      vector3{0.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    CudaPropagatorValidation3, CudaPropConstBFieldMng,
    ::testing::Values(std::make_tuple(-400.f * unit<float>::um,
                                      std::numeric_limits<float>::max(),
                                      vector3{1.f * unit<scalar>::T,
                                              0.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    CudaPropagatorValidation4, CudaPropConstBFieldMng,
    ::testing::Values(std::make_tuple(-600.f * unit<float>::um,
                                      std::numeric_limits<float>::max(),
                                      vector3{1.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    CudaPropagatorValidation5, CudaPropConstBFieldCpy,
    ::testing::Values(std::make_tuple(-100.f * unit<float>::um,
                                      std::numeric_limits<float>::max(),
                                      vector3{0.f * unit<scalar>::T,
                                              0.f * unit<scalar>::T,
                                              2.f * unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    CudaPropagatorValidation6, CudaPropConstBFieldCpy,
    ::testing::Values(std::make_tuple(-400.f * unit<float>::um,
                                      std::numeric_limits<float>::max(),
                                      vector3{0.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    CudaPropagatorValidation7, CudaPropConstBFieldCpy,
    ::testing::Values(std::make_tuple(-400.f * unit<float>::um,
                                      std::numeric_limits<float>::max(),
                                      vector3{1.f * unit<scalar>::T,
                                              0.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    CudaPropagatorValidation8, CudaPropConstBFieldCpy,
    ::testing::Values(std::make_tuple(-600.f * unit<float>::um,
                                      std::numeric_limits<float>::max(),
                                      vector3{1.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T})));

/// This tests the device propagation in an inhomogenepus magnetic field
TEST(CudaPropagatorValidation9, inhomogeneous_bfield_cpy) {
  // VecMem memory resource(s)
  vecmem::host_memory_resource host_mr;
  vecmem::cuda::managed_memory_resource mng_mr;
  vecmem::cuda::device_memory_resource dev_mr;

  vecmem::cuda::copy cuda_cpy;

  // Test configuration
  propagator_test_config cfg{};
  cfg.track_generator.phi_steps(10u).theta_steps(10u);
  cfg.track_generator.p_tot(10.f * unit<scalar>::GeV);
  cfg.track_generator.eta_range(-3.f, 3.f);
  cfg.propagation.navigation.search_window = {3u, 3u};

  // Get the magnetic field
  auto field = create_inhom_field<scalar>();

  // Create the toy geometry with inhomogeneous bfield from file
  auto [det, names] = build_toy_detector<test_algebra>(host_mr);

  auto det_buff = detray::get_buffer(det, dev_mr, cuda_cpy);

  run_propagation_test<bfield::cuda::inhom_bknd_t<scalar>>(
      &mng_mr, det, cfg, detray::get_data(det_buff), std::move(field));
}
