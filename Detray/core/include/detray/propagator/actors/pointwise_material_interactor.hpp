// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/actor.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/material/concepts.hpp"
#include "detray/material/detail/material_accessor.hpp"
#include "detray/material/interaction.hpp"
#include "detray/propagator/actors/parameter_updater.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/utils/geometry_utils.hpp"
#include "detray/utils/logging.hpp"

namespace detray::actor {

template <concepts::algebra algebra_t>
struct pointwise_material_interactor : public base_actor {
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using transform3_type = dtransform3D<algebra_t>;
  using interaction_type = interaction<scalar_type>;
  using bound_param_vector_type = bound_parameters_vector<algebra_t>;
  using bound_matrix_type = bound_matrix<algebra_t>;

  struct state {
    /// Evaluated energy loss
    scalar_type e_loss{0.f};
    /// Evaluated projected scattering angle
    scalar_type projected_scattering_angle{0.f};
    /// Evaluated sigma of qoverp
    scalar_type sigma_qop{0.f};

    bool do_covariance_transport = true;
    bool do_energy_loss = true;
    bool do_multiple_scattering = true;

    DETRAY_HOST_DEVICE
    void reset() {
      e_loss = 0.f;
      projected_scattering_angle = 0.f;
      sigma_qop = 0.f;
    }
  };

  /// Material store visitor
  struct kernel {
    using state = typename pointwise_material_interactor::state;

    template <typename mat_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline bool operator()(
        [[maybe_unused]] const mat_group_t &material_group,
        [[maybe_unused]] const index_t &mat_index, [[maybe_unused]] state &s,
        [[maybe_unused]] const pdg_particle<scalar_type> &ptc,
        [[maybe_unused]] const bound_track_parameters<algebra_t> &bound_params,
        [[maybe_unused]] const scalar_type cos_inc_angle,
        [[maybe_unused]] const scalar_type approach) const {
      using material_t = typename mat_group_t::value_type;

      // Filter material types for which to do "pointwise" interactions
      if constexpr (concepts::surface_material<material_t>) {
        const auto mat = detail::material_accessor::get(
            material_group, mat_index, bound_params.bound_local());

        // Return early in case of zero thickness
        if (mat.thickness() <= std::numeric_limits<scalar_type>::epsilon()) {
          DETRAY_DEBUG_HOST_DEVICE("Zero thickness material: skipping");
          return false;
        }

        DETRAY_VERBOSE_HOST_DEVICE(
            "-> Material path per X0: %f",
            mat.path_segment_in_X0(cos_inc_angle, approach));

        const scalar_type qop = bound_params.qop();

        // Set to exactly zero if the particle is not charged
        if (qop == 0.f) {
          // TODO: Implement interactions for neutral particles here

          DETRAY_DEBUG_HOST_DEVICE(
              "Currently no interactions for neutral particles "
              "implemented: skipping");

          return false;
        }

        const scalar_type path_segment{
            mat.path_segment(cos_inc_angle, approach)};

        detail::relativistic_quantities rq(ptc, qop);

        // Energy Loss
        if (s.do_energy_loss) {
          s.e_loss = interaction_type().compute_energy_loss_bethe_bloch(
              path_segment, mat.get_material(), ptc, rq);
        }

        // @todo: include the radiative loss (Bremsstrahlung)
        if (s.do_energy_loss && s.do_covariance_transport) {
          s.sigma_qop =
              interaction_type().compute_energy_loss_landau_sigma_QOverP(
                  path_segment, mat.get_material(), ptc, rq);
        }

        // Covariance update
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
  DETRAY_HOST_DEVICE inline void operator()(
      state &interactor_state, propagator_state_t &prop_state,
      parameter_transporter_result<algebra_t> &res) const {
    const auto &navigation = prop_state.navigation();

    // Do material interaction when the track is on material surface
    if (!navigation.encountered_sf_material()) {
      return;
    }

    DETRAY_VERBOSE_HOST_DEVICE("Actor: Resolve material effects:");

    interactor_state.reset();

    const auto &stepping = prop_state.stepping();

    const bool success = this->update(
        prop_state.context(), stepping.particle_hypothesis(),
        res.destination_params(), interactor_state,
        static_cast<int>(navigation.direction()), navigation.current_surface());

    res.status = success ? actor::status::e_success : res.status;
  }

  /// @brief Update the bound track parameter
  ///
  /// @param[out] bound_params bound track parameter
  /// @param[out] interactor_state actor state
  /// @param[in]  nav_dir navigation direction
  /// @param[in]  sf the surface
  template <typename context_t, typename surface_t>
  DETRAY_HOST_DEVICE inline bool update(
      const context_t gctx, const pdg_particle<scalar_type> &ptc,
      bound_track_parameters<algebra_t> &bound_params, state &interactor_state,
      const int nav_dir, const surface_t &sf) const {
    // Closest approach of the track to a line surface. Otherwise this is
    // ignored.
    const auto approach{bound_params[e_bound_loc0]};
    const scalar_type cos_inc_angle{
        cos_angle(gctx, sf, bound_params.dir(), bound_params.bound_local())};

    const bool success = sf.template visit_material<kernel>(
        interactor_state, ptc, bound_params, cos_inc_angle, approach);

    if (success) {
      auto &covariance = bound_params.covariance();

      if (interactor_state.do_energy_loss) {
        update_qop(bound_params, ptc, interactor_state.e_loss, nav_dir);

        if (interactor_state.do_covariance_transport) {
          update_qop_variance(covariance, interactor_state.sigma_qop);
        }
      }

      if (interactor_state.do_covariance_transport) {
        update_angle_variance(covariance, bound_params.dir(),
                              interactor_state.projected_scattering_angle);
      }
    }

    assert(!bound_params.is_invalid());

    return success;
  }

  /// @brief Update the q over p of bound track parameter
  ///
  /// @param[out] vector vector of bound track parameter
  /// @param[in]  p momentum of the track
  /// @param[in]  q charge of the track
  /// @param[in]  m mass of the track
  /// @param[in]  e_loss energy loss
  /// @param[in]  sign navigation direction
  DETRAY_HOST_DEVICE
  inline void update_qop(bound_param_vector_type &vector,
                         const pdg_particle<scalar_type> &ptc,
                         const scalar_type e_loss, const int sign) const {
    const scalar_type m = ptc.mass();
    const scalar_type q = ptc.charge();
    const scalar_type p = vector.p(q);

    // Get new Energy
    const scalar_type next_e{
        math::sqrt(m * m + p * p) -
        math::copysign(e_loss, static_cast<scalar_type>(sign))};

    DETRAY_DEBUG_HOST_DEVICE("-> new energy: %f GeV", next_e);

    // Put particle at rest if energy loss is too large
    const scalar_type next_p{(m < next_e) ? math::sqrt(next_e * next_e - m * m)
                                          : 0.f};

    // For neutral particles, qoverp = 1/p
    constexpr auto inv{detail::invalid_value<scalar_type>()};
    const scalar_type next_qop{(q != 0.f) ? q / next_p : 1.f / next_p};
    vector.set_qop((next_p == 0.f) ? inv : next_qop);

    DETRAY_DEBUG_HOST_DEVICE("-> new abs. momentum: %f GeV", next_p);

    DETRAY_DEBUG_HOST_DEVICE("-> Update qop: before: %f e/GeV, after: %f e/GeV",
                             q / p, vector.qop());
  }

  /// @brief Update the variance of q over p of bound track parameter
  ///
  /// @param[out] covariance covariance matrix of bound track parameter
  /// @param[in]  sigma_qop variance of q over p
  /// @param[in]  sign navigation direction
  DETRAY_HOST_DEVICE inline void update_qop_variance(
      bound_matrix_type &covariance, const scalar_type sigma_qop) const {
    const scalar_type variance_qop{sigma_qop * sigma_qop};

    DETRAY_DEBUG_HOST_DEVICE("-> qop variance: %f (e/GeV)^2", variance_qop);

    getter::element(covariance, e_bound_qoverp, e_bound_qoverp) += variance_qop;

    DETRAY_DEBUG_HOST_DEVICE(
        "-> Update qop variance: before: %f (e/GeV)^2, after: %f "
        "(e/GeV)^2",
        getter::element(covariance, e_bound_qoverp, e_bound_qoverp) -
            variance_qop,
        getter::element(covariance, e_bound_qoverp, e_bound_qoverp));
  }

  /// @brief Update the variance of phi and theta of bound track parameter
  ///
  /// @param[out] covariance covariance matrix of bound track parameter
  /// @param[in]  dir direction of track
  /// @param[in]  projected_scattering_angle projected scattering angle
  /// @param[in]  sign navigation direction
  DETRAY_HOST_DEVICE inline void update_angle_variance(
      bound_matrix_type &covariance, const vector3_type &dir,
      const scalar_type projected_scattering_angle) const {
    const scalar_type var_scattering_angle{projected_scattering_angle *
                                           projected_scattering_angle};

    DETRAY_DEBUG_HOST_DEVICE("-> Scattering angle variance: %f rad^2",
                             var_scattering_angle);

    constexpr auto inv{detail::invalid_value<scalar_type>()};

    DETRAY_DEBUG_HOST_DEVICE("-> Update phi variance:");
    DETRAY_DEBUG_HOST_DEVICE(
        "--> before: %f rad^2",
        getter::element(covariance, e_bound_phi, e_bound_phi));

    getter::element(covariance, e_bound_phi, e_bound_phi) +=
        (dir[2] == 1.f) ? inv : var_scattering_angle / (1.f - dir[2] * dir[2]);

    DETRAY_DEBUG_HOST_DEVICE(
        "--> after: %f rad^2",
        getter::element(covariance, e_bound_phi, e_bound_phi));

    DETRAY_DEBUG_HOST_DEVICE("-> Update theta variance:");
    DETRAY_DEBUG_HOST_DEVICE(
        "--> before: %f rad^2",
        getter::element(covariance, e_bound_theta, e_bound_theta));

    getter::element(covariance, e_bound_theta, e_bound_theta) +=
        var_scattering_angle;

    DETRAY_DEBUG_HOST_DEVICE(
        "--> after: %f rad^2",
        getter::element(covariance, e_bound_theta, e_bound_theta));
  }
};

}  // namespace detray::actor
