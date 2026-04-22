// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/pdg_particle.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/definitions/units.hpp"
#include "detray/material/concepts.hpp"
#include "detray/material/detail/material_accessor.hpp"
#include "detray/material/interaction.hpp"
#include "detray/propagator/actors/parameter_updater.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/utils/axis_rotation.hpp"
#include "detray/utils/geometry_utils.hpp"
#include "detray/utils/unit_vectors.hpp"

// Detray test include(s)
#include "detray/test/utils/landau_distribution.hpp"
#include "detray/test/utils/scattering_helper.hpp"

// System include(s).
#include <random>

namespace detray::actor {

template <concepts::algebra algebra_t>
struct random_scatterer : public base_actor {
  using scalar_type = dscalar<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using transform3_type = dtransform3D<algebra_t>;
  using interaction_type = interaction<scalar_type>;

  struct state {
    std::random_device rd{};
    std::mt19937_64 generator{rd()};

    /// most probable energy loss
    scalar_type e_loss_mpv = 0.f;

    /// energy loss sigma
    scalar_type e_loss_sigma = 0.f;

    /// projected scattering angle
    scalar_type projected_scattering_angle = 0.f;

    // Simulation setup
    bool do_energy_loss = true;
    bool do_multiple_scattering = true;

    /// Constructor with seed
    ///
    /// @param sd the seed number
    explicit state(const uint_fast64_t sd = 0u) { generator.seed(sd); }

    void set_seed(const uint_fast64_t sd) { generator.seed(sd); }
  };

  /// Material store visitor
  struct kernel {
    template <typename mat_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline bool operator()(
        [[maybe_unused]] const mat_group_t& material_group,
        [[maybe_unused]] const index_t& mat_index,
        [[maybe_unused]] typename random_scatterer::state& s,
        [[maybe_unused]] const pdg_particle<scalar_type>& ptc,
        [[maybe_unused]] const bound_track_parameters<algebra_t>& bound_params,
        [[maybe_unused]] const scalar_type cos_inc_angle,
        [[maybe_unused]] const scalar_type approach) const {
      using material_t = typename mat_group_t::value_type;

      if constexpr (concepts::surface_material<material_t>) {
        const scalar_type qop = bound_params.qop();

        const auto mat = detail::material_accessor::get(
            material_group, mat_index, bound_params.bound_local());

        // return early in case of zero thickness
        if (mat.thickness() <= std::numeric_limits<scalar_type>::epsilon()) {
          return false;
        }

        const scalar_type path_segment{
            mat.path_segment(cos_inc_angle, approach)};

        const detail::relativistic_quantities rq(ptc, qop);

        // Energy Loss
        if (s.do_energy_loss) {
          s.e_loss_mpv = interaction_type().compute_energy_loss_landau(
              path_segment, mat.get_material(), ptc, rq);

          s.e_loss_sigma = interaction_type().compute_energy_loss_landau_sigma(
              path_segment, mat.get_material(), ptc, rq);
        }

        // Scattering angle
        if (s.do_multiple_scattering) {
          // @todo: use momentum before or after energy loss in
          // backward mode?
          s.projected_scattering_angle =
              interaction_type().compute_multiple_scattering_theta0(
                  mat.path_segment_in_X0(cos_inc_angle, approach), ptc, rq);
        }

        return true;
      } else {
        // For non-pointwise material interactions, do nothing
        return false;
      }
    }
  };

  template <typename propagator_state_t>
  DETRAY_HOST inline void operator()(
      state& simulator_state, propagator_state_t& prop_state,
      parameter_transporter_result<algebra_t>& res) const {
    // @Todo: Make context part of propagation state
    using detector_type = typename propagator_state_t::detector_type;
    using geo_context_type = typename detector_type::geometry_context;

    auto& navigation = prop_state.navigation();

    if (!navigation.encountered_sf_material()) {
      return;
    }

    auto& stepping = prop_state.stepping();
    const auto& ptc = stepping.particle_hypothesis();
    auto& bound_params = res.destination_params();
    const auto sf = navigation.current_surface();
    const scalar_type cos_inc_angle{cos_angle(geo_context_type{}, sf,
                                              bound_params.dir(),
                                              bound_params.bound_local())};

    const bool success = sf.template visit_material<kernel>(
        simulator_state, ptc, bound_params, cos_inc_angle,
        bound_params.bound_local()[0]);

    if (success) {
      // Get the new momentum
      const auto new_mom = attenuate(
          simulator_state.e_loss_mpv, simulator_state.e_loss_sigma, ptc.mass(),
          bound_params.p(ptc.charge()), simulator_state.generator);

      // Update Qop
      bound_params.set_qop(ptc.charge() / new_mom);

      // Get the new direction from random scattering
      const auto new_dir = scatter(bound_params.dir(),
                                   simulator_state.projected_scattering_angle,
                                   simulator_state.generator);

      // Update Phi and Theta
      bound_params.set_phi(vector::phi(new_dir));
      bound_params.set_theta(vector::theta(new_dir));
      // Signal parameter update
      res.status = actor::status::e_success;

      // Flag renavigation of the current candidate
      prop_state.navigation().set_high_trust();
    }
  }

  /// @brief Get the new momentum from the landau distribution
  template <typename generator_t>
  DETRAY_HOST inline scalar_type attenuate(const scalar_type mpv,
                                           const scalar_type sigma,
                                           const scalar_type m0,
                                           const scalar_type p0,
                                           generator_t& generator) const {
    // Get the random energy loss
    // @todo tune the scale parameters (e_loss_mpv and e_loss_sigma)
    const auto e_loss =
        landau_distribution<scalar_type>{}(generator, mpv, sigma);

    // E = sqrt(m^2 + p^2)
    const auto energy = math::sqrt(m0 * m0 + p0 * p0);
    const auto new_energy = energy - e_loss;

    auto p2 = new_energy * new_energy - m0 * m0;

    // To avoid divergence
    if (p2 < 0) {
      p2 = 1.f * unit<scalar_type>::eV * unit<scalar_type>::eV;
    }

    // p = sqrt(E^2 - m^2)
    return math::sqrt(p2);
  }

  /// @brief Scatter the direction with projected scattering angle
  ///
  /// @param dir  input direction
  /// @param projected_scattering_angle  projected scattering angle
  /// @param generator random generator
  /// @returns the new direction from random scattering
  template <typename generator_t>
  DETRAY_HOST inline vector3_type scatter(
      const vector3_type& dir, const scalar_type projected_scattering_angle,
      generator_t& generator) const {
    const auto scattering_angle =
        constant<scalar_type>::sqrt2 * projected_scattering_angle;

    return scattering_helper<algebra_t>()(dir, scattering_angle, generator);
  }
};

}  // namespace detray::actor
