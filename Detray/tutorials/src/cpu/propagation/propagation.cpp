// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/pdg_particle.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/logging.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/track_generators.hpp"

// Example linear algebra plugin: std::array
#include "detray/tutorial/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <sstream>

/// Run the host-propagation
int main() {
  // Toy detector
  using metadata_t = detray::tutorial::toy_metadata;
  using toy_detector_t = detray::detector<metadata_t>;

  using algebra_t = metadata_t::algebra_type;
  using scalar = detray::tutorial::scalar;

  // Navigation
  using navigator_t = detray::caching_navigator<toy_detector_t>;
  // Runge-Kutta-Nystrom stepper (field integration)
  using bfield_t = covfie::field<detray::bfield::const_bknd_t<scalar>>;
  using stepper_t = detray::rk_stepper<bfield_t::view_t, algebra_t>;

  // Actors
  using actor_chain_t = detray::actor_chain<
      detray::actor::pathlimit_aborter<scalar>,
      detray::actor::parameter_updater<
          algebra_t, detray::actor::pointwise_material_interactor<algebra_t>>>;

  // Propagator with empty actor chain
  using propagator_t =
      detray::propagator<stepper_t, navigator_t, actor_chain_t>;

  vecmem::host_memory_resource host_mr;

  std::clog << "Propagation Tutorial\n====================\n\n";
  std::clog << "Building toy detector:\n" << std::endl;

  const auto [det, _] = detray::build_toy_detector<algebra_t>(host_mr);

  typename toy_detector_t::geometry_context gctx{};

  // Create the bfield
  detray::tutorial::vector3 B{0.f, 0.f, 2.f * detray::unit<scalar>::T};
  const bfield_t bfield = detray::create_const_field<scalar>(B);
  const bfield_t::view_t bfield_view(bfield);

  // Build the propagator
  detray::propagation::config prop_cfg{};
  // Particle hypothesis (default: muon)
  constexpr auto ptc = detray::muon<scalar>{};

  propagator_t prop{prop_cfg};

  // Track generation config
  using track_t = detray::free_track_parameters<algebra_t>;
  using track_generator_t = detray::random_track_generator<track_t>;

  track_generator_t::configuration trck_cfg{};
  trck_cfg.n_tracks(1000);
  trck_cfg.eta_range(-3.f, 3.f);
  trck_cfg.pT_range(1.f * detray::unit<scalar>::GeV,
                    100.f * detray::unit<scalar>::GeV);

  std::clog << prop_cfg;
  std::clog << trck_cfg << std::endl;

  // Iterate through uniformly distributed momentum directions
  bool success{true};
  prop_cfg.context = gctx;
  for (auto track : track_generator_t{trck_cfg}) {
    propagator_t::state propagation(track, bfield_view, det, prop_cfg.context);

    propagation.set_particle(detray::update_particle_hypothesis(ptc, track));

    // Prepare actor states
    detray::actor::pathlimit_aborter<scalar>::state aborter_state{
        5.f * detray::unit<scalar>::m};
    detray::actor::parameter_updater_state<algebra_t> updater_state{prop_cfg};
    detray::actor::pointwise_material_interactor<algebra_t>::state
        interactor_state{};

    auto actor_states =
        detray::tie(aborter_state, updater_state, interactor_state);

    // Run the actual propagation
    prop.propagate(propagation, actor_states);
    success &= prop.finished(propagation);
  }

  if (success) {
    DETRAY_INFO_HOST("Successfully propagated " << trck_cfg.n_tracks()
                                                << " tracks!");
  } else {
    DETRAY_ERROR_HOST("Propagation did not complete successfully!");
  }
}
