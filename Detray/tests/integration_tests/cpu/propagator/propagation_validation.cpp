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
#include "detray/propagator/actors.hpp"
#include "detray/propagator/propagation_config.hpp"
#include "detray/tracks/tracks.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/track_generators.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/data_record.hpp"
#include "detray/test/validation/detector_scanner.hpp"
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
// Input tuple: < truth_pT, % max missed surfaces, % max additional surfaces >
class PropagationValidation
    : public ::testing::TestWithParam<
          std::tuple<scalar, unsigned int, float, float>> {};

TEST_P(PropagationValidation, forward_backward) {
  using detector_t = detector<toy_metadata<test_algebra>>;
  using algebra_t = detector_t::algebra_type;
  using bfield_t = bfield::const_field_t<scalar>;
  using track_t = free_track_parameters<algebra_t>;

  using data_record_t = intersection_record<detector_t>;
  using intersection_trace_t = dvector<data_record_t>;
  using generator_t = random_track_generator<track_t>;

  vecmem::host_memory_resource host_mr;

  // Build detector and magnetic field
  toy_det_config<scalar> toy_cfg{};
  toy_cfg.n_brl_layers(4u).n_edc_layers(7u);
  // No material to prevent energy loss during parameter transport
  // (comparing to truth helices)
  toy_cfg.use_material_maps(false).use_homogeneous_material(false);
  const auto [det, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

  // Create b field
  vector3 B{0.f, 0.f, 2.f * unit<scalar>::T};
  const bfield_t hom_bfield = create_const_field<scalar>(B);
  std::optional<bfield_t::view_t> field_view{hom_bfield};

  // Geometry context
  detector_t::geometry_context gctx{};

  // Let the Newton algorithm dynamically choose tol. based on approx. error
  constexpr scalar truth_mask_tol{detray::detail::invalid_value<scalar>()};
  const scalar truth_pT{std::get<0>(GetParam())};

  propagation::config prop_cfg{};
  prop_cfg.navigation.estimate_scattering_noise = false;
  prop_cfg.navigation.search_window = {3u, 3u};

  propagation_validation_config<scalar> test_cfg{};
  test_cfg.propagation = prop_cfg;
  test_cfg.particle = muon<scalar>();
  test_cfg.display_svg = false;  //< faster runtime in the CI
  test_cfg.max_percent_missed = std::get<2>(GetParam());
  test_cfg.max_percent_additional = std::get<3>(GetParam());

  generator_t::configuration trk_gen_cfg{};
  trk_gen_cfg.n_tracks(std::get<1>(GetParam()));
  trk_gen_cfg.p_T(truth_pT);
  trk_gen_cfg.randomize_charge(true);
  // Make sure at least one sensitive is found (otherwise truth trace is empty)
  trk_gen_cfg.eta_range(-4.f, 4.f);

  // Generate the tracks and truth traces for the comparison
  std::vector<track_t> tracks{};
  std::vector<dvector<data_record_t>> truth_traces_fw{};

  for (auto track : generator_t{trk_gen_cfg}) {
    assert(track.qop() != 0.f);
    tracks.push_back(track);

    detail::helix<algebra_t> h{track, B};
    intersection_trace_t intersection_trace = detector_scanner::run<helix_scan>(
        gctx, det, h, truth_mask_tol, track.p(test_cfg.particle.charge()));

    // Only keep the sensitive surfaces
    intersection_trace_t module_trace{};
    for (const data_record_t& rec : intersection_trace) {
      if (rec.intersection.surface().is_sensitive()) {
        module_trace.push_back(rec);
      }
    }
    truth_traces_fw.push_back(std::move(module_trace));
  }

  // Run the propagation and compare the collected data with the truth traces
  const bool success = propagation_validation(det, names, field_view, test_cfg,
                                              tracks, truth_traces_fw);
  ASSERT_TRUE(success);
}

// % of missed surfaces is high likely due to the instability of the
// helix intersections in single precision
INSTANTIATE_TEST_SUITE_P(
    detray_propagator, PropagationValidation,
    ::testing::Values(
        std::make_tuple(100.f * unit<scalar>::GeV, 1000u, 0.2f, 0.1f),
        std::make_tuple(5.f * unit<scalar>::GeV, 1000u, 0.2f, 0.1f),
        std::make_tuple(1.f * unit<scalar>::GeV, 500u, 0.f, 0.f),
        std::make_tuple(0.5f * unit<scalar>::GeV, 500u, 0.2f, 0.1f)));
