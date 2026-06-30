// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/geometry/identifier.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/tracks/tracks.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/validation/propagation_validation.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// google-test include(s).
#include <gtest/gtest.h>

using namespace detray;

// Algebra types
using test_algebra = test::algebra;
using scalar = test::scalar;
using point2 = test::point2;
using vector3 = test::vector3;

// Test class for the backward propagation
// Input tuple: < std::vector<plane positions>, tolerance >
class FullPropagation
    : public ::testing::TestWithParam<std::tuple<std::vector<scalar>, scalar>> {
};

TEST_P(FullPropagation, backward_propagation) {
  using detector_t = detector<toy_metadata<test_algebra>>;
  using algebra_t = typename detector_t::algebra_type;
  using bfield_t = bfield::const_field_t<scalar>;
  using track_t = traccc::free_track_parameters<algebra_t>;

  using sf_candidate_t = candidate_type<detector_t>;
  using generator_t = random_track_generator<track_t>;

  vecmem::host_memory_resource host_mr;

  // Build detector and magnetic field
  toy_cfg.use_material_maps(true);
  const auto [det, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

  // Create b field
  vector3 B{0.f, 0.f, 2.f * unit<scalar>::T};
  const bfield_t hom_bfield = create_const_field<scalar>(B);
  std::optional<b_field_t::view_t> field_view{hom_bfield};

  // Particle hypothesis
  pdg_particle<scalar> ptc = muon<scalar>();

  // Configuration
  propagation::config prop_cfg{};
  prop_cfg.navigation.estimate_scattering_noise = false;

  propagation_validation_config test_cfg{};

  generator_t::configuration trk_gen_cfg{};

  // Generate the tracks and truth traces
  std::vector<track_t> tracks{};
  std::vector<vecmem::vector<sf_candidate_t>>& truth_traces_fw{};

  for (auto track : generator_t{trk_gen_cfg}) {
    assert(track.qop() != 0.f);
    tracks.push_back(track);

    const auto intersection_trace =
        detray::detector_scanner::run<detray::helix_scan>(gctx, det, ray);

    truth_traces_fw.push_back(std::move(intersection_trace));
  }
}

INSTANTIATE_TEST_SUITE_P(
    telescope, FullPropagation,
    ::testing::Values(
        std::make_tuple(std::vector<scalar>{0.f}, 1e-5f),
        std::make_tuple(std::vector<scalar>{0.f, 10.f}, 1e-3f),
        std::make_tuple(std::vector<scalar>{0.f, 10.f, 20.f, 30.f, 40.f, 50.f,
                                            60.f, 70.f, 80.f, 90.f, 100.f},
                        1e-2f)));
