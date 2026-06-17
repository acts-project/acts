// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/material/interaction.hpp"
#include "detray/material/material.hpp"
#include "detray/navigation/policies.hpp"
#include "detray/propagator/base_stepper.hpp"
#include "detray/propagator/codegen/update_rk_transport_jacobian.hpp"
#include "detray/propagator/transport_jacobian.hpp"
#include "detray/tracks/tracks.hpp"

namespace detray {
enum class rk_stepper_flags : std::uint32_t {
  e_allow_covariance_transport = 1,
  e_allow_field_gradient = 2,
};

/// Runge-Kutta-Nystrom 4th order stepper implementation
///
/// @tparam magnetic_field_t the type of magnetic field
/// @tparam track_t the type of track that is being advanced by the stepper
/// @tparam constraint_ the type of constraints on the stepper
template <typename magnetic_field_t, concepts::algebra algebra_t,
          typename constraint_t = unconstrained_step<dscalar<algebra_t>>,
          typename policy_t = stepper_rk_policy<dscalar<algebra_t>>,
          typename inspector_t = stepping::void_inspector,
          std::uint32_t flags_v =
              (static_cast<std::uint32_t>(
                   rk_stepper_flags::e_allow_covariance_transport) |
               static_cast<std::uint32_t>(
                   rk_stepper_flags::e_allow_field_gradient))>
class rk_stepper final
    : public base_stepper<algebra_t, constraint_t, policy_t, inspector_t> {
  using base_type =
      base_stepper<algebra_t, constraint_t, policy_t, inspector_t>;

 public:
  static constexpr bool uses_gradient = static_cast<bool>(
      flags_v &
      static_cast<std::uint32_t>(rk_stepper_flags::e_allow_field_gradient));

  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using point3_type = dpoint3D<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using transform3_type = dtransform3D<algebra_t>;
  using free_track_parameters_type =
      typename base_type::free_track_parameters_type;
  using bound_track_parameters_type =
      typename base_type::bound_track_parameters_type;
  using magnetic_field_type = magnetic_field_t;
  template <std::size_t ROWS, std::size_t COLS>
  using matrix_type = dmatrix<algebra_t, ROWS, COLS>;
  using transport_jacobian_type = std::conditional_t<
      uses_gradient,
      detail::transport_jacobian_matrix_with_gradient<algebra_type>,
      detail::transport_jacobian_matrix_without_gradient<algebra_type>>;

  rk_stepper() = default;

  struct intermediate_state {
    vector3_type b_first{0.f, 0.f, 0.f};
    vector3_type b_middle{0.f, 0.f, 0.f};
    vector3_type b_last{0.f, 0.f, 0.f};
    // t = tangential direction = dr/ds
    darray<vector3_type, 4u> t;
    // q/p
    darray<scalar_type, 4u> qop;
    // dt/ds = d^2r/ds^2 = q/p ( t X B )
    darray<vector3_type, 4u> dtds;
    // d(q/p)/ds
    darray<scalar_type, 4u> dqopds;
  };

  struct state : public base_type::template state<transport_jacobian_type> {
    static constexpr const stepping::id id = stepping::id::e_rk;

    using base_state =
        typename base_type::template state<transport_jacobian_type>;

    DETRAY_HOST_DEVICE
    state(const free_track_parameters_type& t,
          const magnetic_field_t& mag_field)
        : base_state(t), m_magnetic_field(mag_field) {}

    template <typename detector_t>
    DETRAY_HOST_DEVICE state(const bound_track_parameters_type& bound_params,
                             const magnetic_field_t& mag_field,
                             const detector_t& det,
                             const typename detector_t::geometry_context& ctx)
        : base_state(bound_params, det, ctx), m_magnetic_field(mag_field) {}

    /// @returns the B-field view
    DETRAY_HOST_DEVICE
    const magnetic_field_type& field() const { return m_magnetic_field; }

    /// Set the next step size
    DETRAY_HOST_DEVICE
    inline void set_next_step_size(const scalar_type step) {
      assert(math::fabs(step) >= 1e-4f * unit<scalar_type>::mm);
      m_next_step_size = step;
    }

    /// @returns the next step size to be taken on the following step.
    DETRAY_HOST_DEVICE
    inline scalar_type next_step_size() const { return m_next_step_size; }

    /// Evaluate dtds, where t is the unit tangential direction
    DETRAY_HOST_DEVICE
    vector3_type dtds() const {
      // In case there was no step before
      if (this->path_length() == 0.f) [[unlikely]] {
        const point3_type pos = (*this)().pos();

        const auto bvec_tmp = this->field().at(pos[0], pos[1], pos[2]);
        vector3_type bvec;
        bvec[0u] = bvec_tmp[0u];
        bvec[1u] = bvec_tmp[1u];
        bvec[2u] = bvec_tmp[2u];

        DETRAY_DEBUG_HOST(
            "--> dtds: " << (*this)().qop() *
                                vector::cross((*this)().dir(), bvec));

        return (*this)().qop() * vector::cross((*this)().dir(), bvec);
      }

      return m_dtds_3;
    }

    /// Set dtds, where t is the unit tangential direction
    DETRAY_HOST_DEVICE
    void dtds(const vector3_type& v) { m_dtds_3 = v; }

    /// Evaluate d(qop)/ds
    DETRAY_HOST_DEVICE
    scalar_type dqopds(const material<scalar_type>* vol_mat_ptr) const {
      // In case there was no step before
      if (this->path_length() == 0.f) [[unlikely]] {
        return this->dqopds((*this)().qop(), vol_mat_ptr);
      }

      return m_dqopds_3;
    }

    DETRAY_HOST_DEVICE
    scalar_type dqopds(const scalar_type qop,
                       const material<scalar_type>* vol_mat_ptr) const {
      // d(qop)ds is zero for empty space
      if (vol_mat_ptr == nullptr) {
        return 0.f;
      }

      assert(qop != 0.f);
      assert(std::isfinite(qop));
      const scalar_type q = this->particle_hypothesis().charge();
      const scalar_type p = q / qop;
      const scalar_type mass = this->particle_hypothesis().mass();
      const scalar_type E = math::sqrt(p * p + mass * mass);

      // Compute stopping power
      const scalar_type stopping_power =
          interaction<scalar_type>().compute_stopping_power(
              *vol_mat_ptr, this->particle_hypothesis(), {mass, qop, q});

      // Assert that a momentum is a positive value
      assert(p >= 0.f);
      assert(q != 0.f);

      DETRAY_DEBUG_HOST_DEVICE("--> dqopds: %f",
                               qop * qop * qop * E * stopping_power / (q * q));

      // d(qop)ds, which is equal to (qop) * E * (-dE/ds) / p^2
      // or equal to (qop)^3 * E * (-dE/ds) / q^2
      return qop * qop * qop * E * stopping_power / (q * q);
    }

    /// Set d(qop)/ds
    DETRAY_HOST_DEVICE void dqopds(const scalar_type& s) { m_dqopds_3 = s; }

    /// Evaluate d(d(qop)/ds)dqop
    DETRAY_HOST_DEVICE
    scalar_type d2qopdsdqop(const scalar_type qop,
                            const material<scalar_type>* vol_mat_ptr) const {
      if (vol_mat_ptr == nullptr) {
        return 0.f;
      }

      const scalar_type q = this->particle_hypothesis().charge();
      const scalar_type p = q / qop;
      const scalar_type p2 = p * p;

      const auto& mass = this->particle_hypothesis().mass();
      const scalar_type E2 = p2 + mass * mass;

      // Interaction object
      interaction<scalar_type> I;

      // g = dE/ds = -1 * (-dE/ds) = -1 * stopping power
      const detail::relativistic_quantities<scalar_type> rq(mass, qop, q);
      const scalar_type g =
          -1.f * I.compute_stopping_power(*vol_mat_ptr,
                                          this->particle_hypothesis(), rq);

      // dg/d(qop) = -1 * derivation of stopping power
      const scalar_type dgdqop =
          -1.f * I.derive_stopping_power(*vol_mat_ptr,
                                         this->particle_hypothesis(), rq);

      // d(qop)/ds = - qop^3 * E * g / q^2
      const scalar_type dqopds = this->dqopds(qop, vol_mat_ptr);

      // Check Eq 3.12 of
      // (https://iopscience.iop.org/article/10.1088/1748-0221/4/04/P04016/meta)
      assert(E2 != 0.f);
      assert(g != 0.f);
      return dqopds * (1.f / qop * (3.f - p2 / E2) + 1.f / g * dgdqop);
    }

    /// Call the stepping inspector
    template <typename... Args>
    DETRAY_HOST_DEVICE void run_inspector(
        [[maybe_unused]] const stepping::config& cfg,
        [[maybe_unused]] const char* message,
        [[maybe_unused]] const scalar_type dist,
        [[maybe_unused]] Args&&... args) {
      if constexpr (!std::is_same_v<inspector_t, stepping::void_inspector>) {
        this->inspector()(*this, cfg, message, dist,
                          std::forward<Args>(args)...);
      }

      if constexpr (sizeof...(Args) > 0u) {
        DETRAY_DEBUG_HOST("" << message << "\n"
                             << detray::stepping::print_state(
                                    *this, std::forward<Args>(args)...));
      }
    }

   private:
    vector3_type m_dtds_3;
    scalar_type m_dqopds_3;

    /// Next step size after adaptive step size scaling
    scalar_type m_next_step_size{0.f};

    /// Magnetic field view
    const magnetic_field_type& m_magnetic_field;
  };

  /// Take a step, using an adaptive Runge-Kutta algorithm.
  ///
  /// @param dist_to_next The straight line distance to the next surface
  /// @param stepping The state object of a stepper
  /// @param cfg The stepping configuration
  /// @param do_reset whether to reset the RKN step size to "dist to next"
  ///
  /// @return returning the heartbeat, indicating if the stepping is alive
  DETRAY_HOST_DEVICE bool step(
      const scalar_type dist_to_next, state& stepping,
      const stepping::config& cfg, bool do_reset,
      const material<scalar_type>* vol_mat_ptr = nullptr) const {
    DETRAY_DEBUG_HOST("Before: " << stepping());

    // In case of an overlap do nothing
    if (math::fabs(dist_to_next) < 1.f * unit<float>::um) [[unlikely]] {
      DETRAY_VERBOSE_HOST_DEVICE("Zero stepsize...");

      // Don't allow a too small step size on next step
      if (math::fabs(stepping.next_step_size()) < cfg.min_stepsize) {
        stepping.set_next_step_size(math::copysign(
            static_cast<scalar_type>(cfg.min_stepsize), dist_to_next));
      }
      DETRAY_DEBUG_HOST_DEVICE("Setting next stepsize: %f mm",
                               stepping.next_step_size());

      stepping.run_inspector(cfg, "Step skipped (Overlap): ", dist_to_next);

      DETRAY_DEBUG_HOST("After: " << stepping());
      return true;
    }

    if (!(do_reset ||
          math::fabs(stepping.next_step_size()) >= cfg.min_stepsize)) {
      DETRAY_INFO_HOST("Next stepsize: " << stepping.next_step_size() << ", "
                                         << cfg.min_stepsize);
    }

    // Check navigator and actor results
    assert(math::fabs(dist_to_next) != 0.f);
    assert(do_reset ||
           math::fabs(stepping.next_step_size()) >= cfg.min_stepsize);
    assert(!stepping().is_invalid());

    // Get stepper and navigator states
    auto& magnetic_field = stepping.field();

    if (do_reset) {
      stepping.set_step_size(dist_to_next);
    } else if (stepping.next_step_size() > 0) {
      assert(math::fabs(stepping.next_step_size()) >= cfg.min_stepsize);
      stepping.set_step_size(
          math::min(stepping.next_step_size(), dist_to_next));
    } else {
      assert(math::fabs(stepping.next_step_size()) >= cfg.min_stepsize);
      stepping.set_step_size(
          math::max(stepping.next_step_size(), dist_to_next));
    }

    DETRAY_VERBOSE_HOST_DEVICE("Distance to nex: %f mm", dist_to_next);
    DETRAY_VERBOSE_HOST_DEVICE("Initial stepsize: %f mm", stepping.step_size());

    // Don't allow too small stepsizes, unless the navigation needs it
    const scalar_type min_stepsize{
        math::min(math::fabs(stepping.step_size()),
                  static_cast<scalar_type>(cfg.min_stepsize))};
    const point3_type pos = stepping().pos();

    intermediate_state sd{};

    // First Runge-Kutta point
    auto bvec = magnetic_field.at(pos[0], pos[1], pos[2]);
    assert(math::isfinite(bvec[0]));
    assert(math::isfinite(bvec[1]));
    assert(math::isfinite(bvec[2]));
    sd.b_first[0] = bvec[0];
    sd.b_first[1] = bvec[1];
    sd.b_first[2] = bvec[2];

    DETRAY_DEBUG_HOST_DEVICE("First stage:");
    DETRAY_DEBUG_HOST_DEVICE("-> B-field: [%f, %f, %f]", bvec[0], bvec[1],
                             bvec[2]);

    // qop should be recalculated at every point
    // Reference: Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
    detray::tie(sd.dqopds[0u], sd.qop[0u]) =
        evaluate_dqopds(stepping, 0u, 0.f, 0.f, vol_mat_ptr, cfg);
    detray::tie(sd.dtds[0u], sd.t[0u]) = evaluate_dtds(
        stepping, sd.b_first, 0u, 0.f, vector3_type{0.f, 0.f, 0.f}, sd.qop[0u]);

    /// RKN step trial and error estimation
    const auto estimate_error = [&](const scalar_type& h) {
      assert(h != 0);
      // State the square and half of the step size
      const scalar_type h2{h * h};
      const scalar_type half_h{h * 0.5f};

      // Second Runge-Kutta point
      // qop should be recalculated at every point
      // Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
      const point3_type pos1 =
          pos + half_h * sd.t[0u] + h2 * 0.125f * sd.dtds[0u];
      bvec = magnetic_field.at(pos1[0], pos1[1], pos1[2]);
      assert(math::isfinite(bvec[0]));
      assert(math::isfinite(bvec[1]));
      assert(math::isfinite(bvec[2]));
      sd.b_middle[0] = bvec[0];
      sd.b_middle[1] = bvec[1];
      sd.b_middle[2] = bvec[2];
      DETRAY_DEBUG_HOST_DEVICE("Second stage:");
      DETRAY_DEBUG_HOST_DEVICE("-> B-field: [%f, %f, %f]", bvec[0], bvec[1],
                               bvec[2]);

      detray::tie(sd.dqopds[1u], sd.qop[1u]) = evaluate_dqopds(
          stepping, 1u, half_h, sd.dqopds[0u], vol_mat_ptr, cfg);
      detray::tie(sd.dtds[1u], sd.t[1u]) = evaluate_dtds(
          stepping, sd.b_middle, 1u, half_h, sd.dtds[0u], sd.qop[1u]);

      // Third Runge-Kutta point
      // qop should be recalculated at every point
      // Reference: Eq (84) of
      // https://doi.org/10.1016/0029-554X(81)90063-X
      detray::tie(sd.dqopds[2u], sd.qop[2u]) = evaluate_dqopds(
          stepping, 2u, half_h, sd.dqopds[1u], vol_mat_ptr, cfg);
      detray::tie(sd.dtds[2u], sd.t[2u]) = evaluate_dtds(
          stepping, sd.b_middle, 2u, half_h, sd.dtds[1u], sd.qop[2u]);

      // Last Runge-Kutta point
      // qop should be recalculated at every point
      // Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
      const point3_type pos2 = pos + h * sd.t[0u] + h2 * 0.5f * sd.dtds[2u];
      bvec = magnetic_field.at(pos2[0], pos2[1], pos2[2]);
      assert(math::isfinite(bvec[0]));
      assert(math::isfinite(bvec[1]));
      assert(math::isfinite(bvec[2]));
      sd.b_last[0] = bvec[0];
      sd.b_last[1] = bvec[1];
      sd.b_last[2] = bvec[2];

      DETRAY_DEBUG_HOST_DEVICE("Third stage:");
      DETRAY_DEBUG_HOST_DEVICE("-> B-field: [%f, %f, %f]", bvec[0], bvec[1],
                               bvec[2]);

      detray::tie(sd.dqopds[3u], sd.qop[3u]) =
          evaluate_dqopds(stepping, 3u, h, sd.dqopds[2u], vol_mat_ptr, cfg);
      detray::tie(sd.dtds[3u], sd.t[3u]) =
          evaluate_dtds(stepping, sd.b_last, 3u, h, sd.dtds[2u], sd.qop[3u]);

      // Compute and check the local integration error estimate
      // @Todo
      constexpr auto one_sixth{static_cast<scalar_type>(1. / 6.)};
      const vector3_type err_vec =
          one_sixth * h2 *
          (sd.dtds[0u] - sd.dtds[1u] - sd.dtds[2u] + sd.dtds[3u]);

      DETRAY_DEBUG_HOST("-> Integration error vector: " << err_vec);

      return vector::norm(err_vec);
    };

    /// Calculate the scale factor for the stepsize adjustment using the
    /// error estimate @param err
    const auto step_size_scaling =
        [&cfg](const scalar_type& err) -> scalar_type {
      assert(err != 0.f);
      return static_cast<scalar_type>(
          math::min(math::max(math::sqrt(math::sqrt(cfg.rk_error_tol / err)),
                              static_cast<scalar_type>(0.25)),
                    static_cast<scalar_type>(4.)));
    };

    scalar_type error{1e20f};

    // If the estimated error is larger than the tolerance with an
    // additional margin, reduce the step size and try again
    const auto n_trials{cfg.max_rk_updates};
    for (unsigned int i = 0u; i < n_trials; i++) {
      stepping.count_trials();

      error = math::max(estimate_error(stepping.step_size()),
                        static_cast<scalar_type>(1e-20));

      DETRAY_DEBUG_HOST_DEVICE("-> Integration error: %f", error);
      assert(math::isfinite(error));

      // Error is small enough
      // ---> break and advance track
      if (error <= 4.f * cfg.rk_error_tol) {
        break;
      }
      // Error estimate is too big
      // ---> Make step size smaller and estimate error again
      else {
        stepping.set_step_size(stepping.step_size() * step_size_scaling(error));

        // Run inspection while the stepsize is getting adjusted
        stepping.run_inspector(cfg, "Adjust stepsize: ", dist_to_next, i,
                               step_size_scaling(error));
      }
    }

    stepping.dtds(sd.dtds[3u]);
    stepping.dqopds(sd.dqopds[3u]);

    // Check constraints
    if (const scalar_type max_step =
            stepping.constraints().template size<>(stepping.direction());
        math::fabs(stepping.step_size()) > math::fabs(max_step)) {
      // Run inspection before step size is cut
      stepping.run_inspector(cfg, "Before constraint: ", dist_to_next);

      stepping.set_step_size(max_step);
    }

    // Adjust the min step size
    if (math::fabs(stepping.step_size()) < min_stepsize) {
      stepping.set_step_size(
          math::copysign(min_stepsize, stepping.step_size()));
    }

    DETRAY_VERBOSE_HOST_DEVICE("Take step: %f mm", stepping.step_size());

    // The step size estimation for the next step
    stepping.set_next_step_size(stepping.step_size() *
                                step_size_scaling(error));

    // Don't allow a too small step size
    if (math::fabs(stepping.next_step_size()) < cfg.min_stepsize) {
      stepping.set_next_step_size(
          math::copysign(static_cast<scalar_type>(cfg.min_stepsize),
                         stepping.next_step_size()));
    }

    DETRAY_DEBUG_HOST_DEVICE("Estimated next stepsize: %f mm",
                             stepping.next_step_size());

    assert(math::fabs(stepping.next_step_size()) >= cfg.min_stepsize);
    assert(math::fabs(stepping.step_size()) >= cfg.min_stepsize);

    // Advance track state
    advance_track(stepping, sd, vol_mat_ptr);
    assert(!stepping().is_invalid());

    // Advance jacobian transport
    if constexpr ((flags_v & static_cast<std::uint32_t>(
                                 detray::rk_stepper_flags::
                                     e_allow_covariance_transport)) != 0u) {
      if (cfg.do_covariance_transport) {
        advance_jacobian(stepping, cfg, sd, vol_mat_ptr);
      }
    } else {
      assert(!cfg.do_covariance_transport);
    }

    // Run final inspection
    stepping.run_inspector(cfg, "Step complete: ", dist_to_next);

    DETRAY_DEBUG_HOST("After: " << stepping());
    return true;
  }

  /// evaluate dqopds for a given step size and material
  DETRAY_HOST_DEVICE
  detray::pair<scalar_type, scalar_type> evaluate_dqopds(
      detray::rk_stepper<magnetic_field_t, algebra_t, constraint_t, policy_t,
                         inspector_t, flags_v>::state& stepping,
      const std::size_t i, const scalar_type h, const scalar_type dqopds_prev,
      const material<scalar_type>* vol_mat_ptr,
      const detray::stepping::config& cfg) const {
    const auto& track = stepping();

    if (vol_mat_ptr == nullptr) {
      const scalar_type qop = track.qop();
      DETRAY_DEBUG_HOST_DEVICE("-> qop: %f", qop);
      DETRAY_DEBUG_HOST_DEVICE("-> dqopds: 0");

      return detray::make_pair(scalar_type(0.f), qop);
    } else if (cfg.use_mean_loss && i != 0u) {
      // qop_n is calculated recursively like the direction of
      // evaluate_dtds.
      //
      // https://doi.org/10.1016/0029-554X(81)90063-X says:
      // "For y  we  have  similar  formulae  as  for x, for y' and
      // \lambda similar  formulae as for  x'"
      const scalar_type qop = track.qop() + h * dqopds_prev;
      DETRAY_DEBUG_HOST_DEVICE("-> qop: %f", qop);

      return detray::make_pair(stepping.dqopds(qop, vol_mat_ptr), qop);
    } else {
      const scalar_type qop = track.qop();
      DETRAY_DEBUG_HOST_DEVICE("-> qop: %f", qop);

      return detray::make_pair(stepping.dqopds(qop, vol_mat_ptr), qop);
    }
  }

  /// Evaluate dtds for Runge-Kutta stepping
  DETRAY_HOST_DEVICE
  detray::pair<vector3_type, vector3_type> evaluate_dtds(
      detray::rk_stepper<magnetic_field_t, algebra_t, constraint_t, policy_t,
                         inspector_t, flags_v>::state& stepping,
      const vector3_type& b_field, const std::size_t i, const scalar_type h,
      const vector3_type& dtds_prev, const scalar_type qop) const {
    auto& track = stepping();
    const auto dir = track.dir();

    assert(std::isfinite(qop));

    // Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
    vector3_type t{(i == 0u || h == 0.f)
                       ? dir
                       : static_cast<vector3_type>(dir + h * dtds_prev)};

    // dtds = qop * (t X B) from Lorentz force
    DETRAY_DEBUG_HOST(
        "-> evaluate dtds: " << vector3_type{qop * vector::cross(t, b_field)});

    return detray::make_pair(vector3_type{qop * vector::cross(t, b_field)}, t);
  }

  /// Update the track state by Runge-Kutta-Nystrom integration.
  DETRAY_HOST_DEVICE
  void advance_track(
      detray::rk_stepper<magnetic_field_t, algebra_t, constraint_t, policy_t,
                         inspector_t, flags_v>::state& stepping,
      const intermediate_state& sd,
      const material<scalar_type>* vol_mat_ptr) const {
    const scalar_type h{stepping.step_size()};
    const scalar_type h_6{h * static_cast<scalar_type>(1. / 6.)};
    auto& track = stepping();
    auto pos = track.pos();
    auto dir = track.dir();

    // Update the track parameters according to the equations of motion
    // Reference: Eq (82) of https://doi.org/10.1016/0029-554X(81)90063-X
    pos = pos + h * (sd.t[0u] + h_6 * (sd.dtds[0] + sd.dtds[1] + sd.dtds[2]));
    track.set_pos(pos);

    // Reference: Eq (82) of https://doi.org/10.1016/0029-554X(81)90063-X
    dir =
        dir + h_6 * (sd.dtds[0] + 2.f * (sd.dtds[1] + sd.dtds[2]) + sd.dtds[3]);
    dir = vector::normalize(dir);
    track.set_dir(dir);

    auto qop = track.qop();
    if (vol_mat_ptr != nullptr) {
      // Reference: Eq (82) of
      // https://doi.org/10.1016/0029-554X(81)90063-X
      qop = qop + h_6 * (sd.dqopds[0u] + 2.f * (sd.dqopds[1u] + sd.dqopds[2u]) +
                         sd.dqopds[3u]);
    }
    track.set_qop(qop);

    // Update path lengths
    stepping.update_path_lengths(h);
  }

  /// Update the jacobian transport from free propagation
  DETRAY_HOST_DEVICE
  void advance_jacobian(
      detray::rk_stepper<magnetic_field_t, algebra_t, constraint_t, policy_t,
                         inspector_t, flags_v>::state& stepping,
      const stepping::config& cfg, const intermediate_state& sd,
      const material<scalar_type>* vol_mat_ptr) const {
    DETRAY_VERBOSE_HOST_DEVICE("-> Advance Jacobian...");

    /// The calculations are based on ATL-SOFT-PUB-2009-002. The update of
    /// the Jacobian matrix is requires only the calculation of eq. 17
    /// and 18. Since the terms of eq. 18 are currently 0, this matrix is
    /// not needed in the calculation. The matrix A from eq. 17 consists out
    /// of 3 different parts. The first one is given by the upper left 3x3
    /// matrix that are calculated by the derivatives dF/dT (called dFdT)
    /// and dG/dT (calls dGdT). The second is given by the top 3 lines of
    /// the rightmost column. This is calculated by dFdqop and dGdqop. The
    /// remaining non-zero term is calculated directly. The naming of the
    /// variables is explained in eq. 11 and are directly related to the
    /// initial problem in eq. 7. The evaluation is based by propagating the
    /// parameters T and lambda as given in eq. 16 and evaluating the
    /// derivations for matrix A.
    /// @note The translation for u_{n+1} in eq. 7 is in this case a
    /// 3-dimensional vector without a dependency of Lambda or lambda
    /// neither in u_n nor in u_n'. The second and fourth eq. in eq. 14 have
    /// the constant offset matrices h * Id and Id respectively. This
    /// involves that the constant offset does not exist for rectangular
    /// matrix dGdu' (due to the missing Lambda part) and only exists for
    /// dFdu' in dlambda/dlambda.

    const scalar_type h{stepping.step_size()};
    const scalar_type h2{h * h};
    const scalar_type half_h{h * 0.5f};
    const scalar_type h_6{h * (1.f / 6.f)};

    /*---------------------------------------------------------------------------
     *  dk_n/dt1
     *    = qop_n * (dt_n/dt1 X B_n)
     *      + qop_n * ( t_n X dB_n/dt1 ),
     *  where dB_n/dt1 == dB_n/dr_n * dr_n/dt1.
     *
     *  The second term is non-zero only for inhomogeneous magnetic fields
     *
     *  Note that [ t_n = t1 + h * d(t_{n-1})/ds) ] as indicated by
     *  Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X

     *  [ Table for dt_n/dt1 ]
     *  dt1/dt1 = I
     *  dt2/dt1 = d( t1 + h/2 * dt1/ds ) / dt1 = I + h/2 * dk1/dt1
     *  dt3/dt1 = d( t1 + h/2 * dt2/ds ) / dt1 = I + h/2 * dk2/dt1
     *  dt4/dt1 = d( t1 + h * dt3/ds ) / dt1 = I + h * dk3/dt1
     *
     *  [ Table for dr_n/dt1 ]
     *  dr1/dt1 = 0
     *  dr2/dt1 = d(r1 + h/2 * t1 + h^2/8 dt1/ds)/dt1 = h/2 * I + h^2/8
    dk1/dt1
     *  dr3/dt1 = d(r1 + h/2 * t1 + h^2/8 dt1/ds)/dt1 = h/2 * I + h^2/8
    dk1/dt1
     *  dr4/dt1 = d(r1 + h * t1 + h^2/2 dt3/ds)/dt1 = h * I + h^2/2 dk3/dt1
     *
     *  Note that
     *  d/dr [ F(T) X B ]  = dF(T)/dr (X) B, where (X) means the column wise
     *  cross product
    ---------------------------------------------------------------------------*/

    /*---------------------------------------------------------------------------
     *  dk_n/dqop_1
     *    = dqop_n/dqop1 * ( t_n X B_n )
     *      + qop_n * ( dt_n/dqop1 X B_n )
     *      + qop_n * ( t_n X dB_n/dqop1 ),
     *  where dB_n/dqop1 = dB_n/dr_n * dr_n/dqop1
     *
     *  Note that [ qop_n = qop1 + h * dqop_{n-1}/ds) ] as indicated by
     *  Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
     *
     *  [ Table for dqop_n/dqop1 ]
     *  dqop1/dqop1 = 1
     *  dqop2/dqop1 = 1 + h/2 * d(dqop1/ds)/dqop1
     *  dqop3/dqop1 = 1 + h/2 * d(dqop2/ds)/dqop1
     *  dqop4/dqop1 = 1 + h * d(dqop3/ds)/dqop1
     *
     *  [ Table for dt_n/dqop1 ]
     *  dt1/dqop1 = 0
     *  dt2/dqop1 = d(t1 + h/2 dt1/ds)/dqop1 = h/2 * dk1/dqop1
     *  dt3/dqop1 = d(t1 + h/2 dt2/ds)/dqop1 = h/2 * dk2/dqop1
     *  dt4/dqop1 = d(t1 + h dt3/ds)/dqop1 = h * dk3/dqop1
     *
     *  [ Table for dr_n/dqop1 ]
     *  dr1/dqop1 = 0
     *  dr2/dqop1 = d(r1 + h/2 * t1 + h^2/8 dt1/ds)/dqop1 = h^2/8 *
    dk1/dqop1
     *  dr3/dqop1 = d(r1 + h/2 * t1 + h^2/8 dt1/ds)/dqop1 = h^2/8 *
    dk1/dqop1
     *  dr4/dqop1 = d(r1 + h * t1 + h^2/2 dt3/ds)/dqop1 = h^2/2 dk3/dqop1
    ---------------------------------------------------------------------------*/

    /*---------------------------------------------------------------------------
     *  dk_n/dr1
     *    = qop_n * ( dt_n/dr1 X B_n )
     *      + qop_n * ( t_n X dB_n/dr1 ),
     *  where dB_n/dr1 = dB_n/dr_n * dr_n/dr1
     *
     *  [ Table for dt_n/dr1 ]
     *  dt1/dr1 = 0
     *  dt2/dr1 = d(t1 + h/2 * dt1/ds)/dr1 = h/2 * dk1/dr1
     *  dt2/dr1 = d(t1 + h/2 * dt2/ds)/dr1 = h/2 * dk2/dr1
     *  dt3/dr1 = d(t1 + h * dt3/ds)/dr1 = h * dk3/dr1
     *
     *  [ Table for dr_n/dr1 ]
     *  dr1/dr1 = I
     *  dr2/dr1 = (r1 + h/2 * t1 + h^2/8 dt1/ds ) / dr1 = I + h^2/8 dk1/dr1
     *  dr3/dr1 = (r1 + h/2 * t1 + h^2/8 dt1/ds ) / dr1 = I + h^2/8 dk1/dr1
     *  dr4/dr1 = (r1 + h * t1 + h^2/2 dt3/ds ) / dr1 = I + h^2/2 dk3/dr1
    ---------------------------------------------------------------------------*/

    /*---------------------------------------------------------------------------
     *  d(dqop_n/ds)/dqop1
     *
     *  Useful equation:
     *  dqop/ds = qop^3 * E * (-dE/ds) / q^2 = - qop^3 * E g / q^2

     *  [ Table for d(dqop_n/ds)/dqop1 ]
     *  d(dqop1/ds)/dqop1 = dqop1/ds * (1/qop * (3 - p^2/E^2) + 1/g1 *
    dg1dqop1
     *  d(dqop2/ds)/dqop1 = d(dqop2/ds)/dqop2 * (1 + h/2 *
    d(dqop1/ds)/dqop1)
     *  d(dqop3/ds)/dqop1 = d(dqop3/ds)/dqop3 * (1 + h/2 *
    d(dqop2/ds)/dqop1)
     *  d(dqop4/ds)/dqop1 = d(dqop4/ds)/dqop4 * (1 + h * d(dqop3/ds)/dqop1)
    ---------------------------------------------------------------------------*/

    scalar_type dqopqop;
    vector3_type dFdqop;
    vector3_type dGdqop;

    {
      darray<scalar_type, 4u> dqopn_dqop{1.f, 1.f, 1.f, 1.f};

      if (!cfg.use_eloss_gradient) {
        dqopqop = 1.f;
      } else {
        // Pre-calculate dqop_n/dqop1
        const scalar_type d2qop1dsdqop1 =
            stepping.d2qopdsdqop(sd.qop[0u], vol_mat_ptr);

        dqopn_dqop[0u] = 1.f;
        dqopn_dqop[1u] = 1.f + half_h * d2qop1dsdqop1;

        const scalar_type d2qop2dsdqop1 =
            stepping.d2qopdsdqop(sd.qop[1u], vol_mat_ptr) * dqopn_dqop[1u];
        dqopn_dqop[2u] = 1.f + half_h * d2qop2dsdqop1;

        const scalar_type d2qop3dsdqop1 =
            stepping.d2qopdsdqop(sd.qop[2u], vol_mat_ptr) * dqopn_dqop[2u];
        dqopn_dqop[3u] = 1.f + h * d2qop3dsdqop1;

        const scalar_type d2qop4dsdqop1 =
            stepping.d2qopdsdqop(sd.qop[3u], vol_mat_ptr) * dqopn_dqop[3u];

        /*-----------------------------------------------------------------
         * Calculate the first terms of d(dqop_n/ds)/dqop1
        -------------------------------------------------------------------*/

        dqopqop =
            1.f + h_6 * (d2qop1dsdqop1 + 2.f * (d2qop2dsdqop1 + d2qop3dsdqop1) +
                         d2qop4dsdqop1);
      }

      /*-----------------------------------------------------------------
       * Calculate the first and second terms of dk_n/dqop1
      -------------------------------------------------------------------*/

      {
        darray<vector3_type, 4u> dkndqop;

        // dk1/dqop1
        dkndqop[0u] = dqopn_dqop[0u] * vector::cross(sd.t[0u], sd.b_first);

        // dk2/dqop1
        dkndqop[1u] =
            dqopn_dqop[1u] * vector::cross(sd.t[1u], sd.b_middle) +
            sd.qop[1u] * half_h * vector::cross(dkndqop[0u], sd.b_middle);

        // dk3/dqop1
        dkndqop[2u] =
            dqopn_dqop[2u] * vector::cross(sd.t[2u], sd.b_middle) +
            sd.qop[2u] * half_h * vector::cross(dkndqop[1u], sd.b_middle);

        // dk4/dqop1
        dkndqop[3u] = dqopn_dqop[3u] * vector::cross(sd.t[3u], sd.b_last) +
                      sd.qop[3u] * h * vector::cross(dkndqop[2u], sd.b_last);

        // Set dF/dqop1 and dG/dqop1
        dFdqop = h * h_6 * (dkndqop[0u] + dkndqop[1u] + dkndqop[2u]);
        dGdqop = h_6 * (dkndqop[0u] + 2.f * (dkndqop[1u] + dkndqop[2u]) +
                        dkndqop[3u]);
      }
    }

    /*-----------------------------------------------------------------
     * Calculate the first terms of dk_n/dt1
    -------------------------------------------------------------------*/
    // Set dF/dt1 and dG/dt1
    auto dFdt = matrix::identity<matrix_type<3, 3>>();
    auto dGdt = matrix::identity<matrix_type<3, 3>>();

    {
      const auto I33 = matrix::identity<matrix_type<3, 3>>();
      darray<matrix_type<3u, 3u>, 4u> dkndt{I33, I33, I33, I33};

      // dk1/dt1
      dkndt[0u] = sd.qop[0u] * matrix::column_wise_cross(dkndt[0u], sd.b_first);

      // dk2/dt1
      dkndt[1u] = dkndt[1u] + half_h * dkndt[0u];
      dkndt[1u] =
          sd.qop[1u] * matrix::column_wise_cross(dkndt[1u], sd.b_middle);

      // dk3/dt1
      dkndt[2u] = dkndt[2u] + half_h * dkndt[1u];
      dkndt[2u] =
          sd.qop[2u] * matrix::column_wise_cross(dkndt[2u], sd.b_middle);

      // dk4/dt1
      dkndt[3u] = dkndt[3u] + h * dkndt[2u];
      dkndt[3u] = sd.qop[3u] * matrix::column_wise_cross(dkndt[3u], sd.b_last);

      dFdt = dFdt + h_6 * (dkndt[0u] + dkndt[1u] + dkndt[2u]);
      dFdt = h * dFdt;
      dGdt =
          dGdt + h_6 * (dkndt[0u] + 2.f * (dkndt[1u] + dkndt[2u]) + dkndt[3u]);
    }

    // Calculate dkndr in case of considering B field gradient
    auto dFdr = matrix::identity<matrix_type<3, 3>>();
    auto dGdr = matrix::zero<matrix_type<3, 3>>();

    if constexpr ((flags_v & static_cast<std::uint32_t>(
                                 rk_stepper_flags::e_allow_field_gradient)) !=
                  0u) {
      if (cfg.use_field_gradient) {
        darray<matrix_type<3u, 3u>, 4u> dkndr;
        auto& track = stepping();

        // Positions and field gradients at initial, middle and final
        // points of the fourth order RKN
        vector3_type r_ini = track.pos();
        vector3_type r_mid =
            r_ini + half_h * sd.t[0u] + h2 * 0.125f * sd.dtds[0u];
        vector3_type r_fin = r_ini + h * sd.t[0u] + h2 * 0.5f * sd.dtds[2u];

        matrix_type<3, 3> dBdr_ini = evaluate_field_gradient(stepping, r_ini);
        matrix_type<3, 3> dBdr_mid = evaluate_field_gradient(stepping, r_mid);
        matrix_type<3, 3> dBdr_fin = evaluate_field_gradient(stepping, r_fin);

        /*-----------------------------------------------------------------
         * Calculate all terms of dk_n/dr1
        -------------------------------------------------------------------*/
        const auto I33 = matrix::identity<matrix_type<3, 3>>();

        // dk1/dr1
        dkndr[0u] = -sd.qop[0u] * matrix::column_wise_cross(dBdr_ini, sd.t[0u]);

        const auto dkndr0_tmp = (I33 + h2 * 0.125f * dkndr[0u]);

        // dk2/dr1
        dkndr[1u] = sd.qop[1u] *
                    matrix::column_wise_cross(half_h * dkndr[0u], sd.b_middle);
        dkndr[1u] =
            dkndr[1u] - sd.qop[1u] * matrix::column_wise_cross(
                                         dBdr_mid * dkndr0_tmp, sd.t[1u]);

        // dk3/dr1
        dkndr[2u] = sd.qop[2u] *
                    matrix::column_wise_cross(half_h * dkndr[1u], sd.b_middle);
        dkndr[2u] =
            dkndr[2u] - sd.qop[2u] * matrix::column_wise_cross(
                                         dBdr_mid * dkndr0_tmp, sd.t[2u]);

        // dk4/dr1
        dkndr[3u] =
            sd.qop[3u] * matrix::column_wise_cross(h * dkndr[2u], sd.b_last);
        dkndr[3u] = dkndr[3u] -
                    sd.qop[3u] *
                        matrix::column_wise_cross(
                            dBdr_fin * (I33 + h2 * 0.5f * dkndr[2u]), sd.t[3u]);

        // Set dF/dr1 and dG/dr1
        dFdr = dFdr + h * h_6 * (dkndr[0u] + dkndr[1u] + dkndr[2u]);
        dGdr = h_6 * (dkndr[0u] + 2.f * (dkndr[1u] + dkndr[2u]) + dkndr[3u]);
      }
    } else {
      assert(!cfg.use_field_gradient);
    }

    const auto old_jacobian = stepping.internal_transport_jacobian();

    if constexpr ((flags_v & static_cast<std::uint32_t>(
                                 rk_stepper_flags::e_allow_field_gradient)) !=
                  0u) {
      detail::update_transport_jacobian_with_gradient_impl(
          old_jacobian, dFdt, dGdt, dFdr, dGdr, dFdqop, dGdqop, dqopqop,
          stepping.internal_transport_jacobian());
    } else {
      assert(!cfg.use_field_gradient);

      detail::update_transport_jacobian_without_gradient_impl(
          old_jacobian, dFdt, dGdt, dFdqop, dGdqop, dqopqop,
          stepping.internal_transport_jacobian());
    }
  }

  DETRAY_HOST_DEVICE
  matrix_type<3, 3> evaluate_field_gradient(
      detray::rk_stepper<magnetic_field_t, algebra_t, constraint_t, policy_t,
                         inspector_t, flags_v>::state& stepping,
      const point3_type& pos) const {
    auto dBdr = matrix::zero<matrix_type<3, 3>>();

    constexpr auto delta{1e-1f * unit<scalar_type>::mm};

    for (unsigned int i = 0; i < 3; i++) {
      point3_type dpos1 = pos;
      dpos1[i] += delta;
      const auto bvec1_tmp = stepping.field().at(dpos1[0], dpos1[1], dpos1[2]);
      vector3_type bvec1;
      bvec1[0u] = bvec1_tmp[0u];
      bvec1[1u] = bvec1_tmp[1u];
      bvec1[2u] = bvec1_tmp[2u];

      point3_type dpos2 = pos;
      dpos2[i] -= delta;
      const auto bvec2_tmp = stepping.field().at(dpos2[0], dpos2[1], dpos2[2]);
      vector3_type bvec2;
      bvec2[0u] = bvec2_tmp[0u];
      bvec2[1u] = bvec2_tmp[1u];
      bvec2[2u] = bvec2_tmp[2u];

      assert(delta != 0.f);
      const vector3_type gradient = (bvec1 - bvec2) * (1.f / (2.f * delta));

      getter::element(dBdr, 0u, i) = gradient[0u];
      getter::element(dBdr, 1u, i) = gradient[1u];
      getter::element(dBdr, 2u, i) = gradient[2u];
    }

    return dBdr;
  }
};

}  // namespace detray
