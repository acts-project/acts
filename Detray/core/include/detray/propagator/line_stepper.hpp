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
#include "detray/definitions/math.hpp"
#include "detray/material/material.hpp"
#include "detray/navigation/policies.hpp"
#include "detray/propagator/base_stepper.hpp"
#include "detray/utils/logging.hpp"

namespace detray {

/// Straight line stepper implementation
template <concepts::algebra algebra_t,
          typename constraint_t = unconstrained_step<dscalar<algebra_t>>,
          typename policy_t = stepper_default_policy<dscalar<algebra_t>>,
          typename inspector_t = stepping::void_inspector>
class line_stepper final
    : public base_stepper<algebra_t, constraint_t, policy_t, inspector_t> {
  using base_type =
      base_stepper<algebra_t, constraint_t, policy_t, inspector_t>;

 public:
  static constexpr bool uses_gradient = true;

  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using transform3_type = dtransform3D<algebra_t>;
  using free_track_parameters_type =
      typename base_type::free_track_parameters_type;
  using bound_track_parameters_type =
      typename base_type::bound_track_parameters_type;
  template <std::size_t ROWS, std::size_t COLS>
  using matrix_type = dmatrix<algebra_t, ROWS, COLS>;

  struct state : public base_type::template state<> {
    static constexpr const stepping::id id = stepping::id::e_linear;

    using base_state = typename base_type::template state<>;

    /// Import base class constructors
    using base_state::base_state;

    /// Update the track state in a straight line.
    DETRAY_HOST_DEVICE
    inline void advance_track() {
      auto& track = (*this)();
      track.set_pos(track.pos() + track.dir() * this->step_size());

      this->update_path_lengths(this->step_size());
    }

    DETRAY_HOST_DEVICE
    inline void advance_jacobian() {
      DETRAY_VERBOSE_HOST_DEVICE("Advance Jacobian");

      // The step transport matrix in global coordinates
      free_matrix<algebra_t> D = matrix::identity<free_matrix<algebra_t>>();

      // d(x,y,z)/d(n_x,n_y,n_z)
      matrix_type<3, 3> dxdn =
          this->step_size() * matrix::identity<matrix_type<3, 3>>();
      getter::set_block(D, dxdn, e_free_pos0, e_free_dir0);

      /// NOTE: Let's skip the element for d(time)/d(qoverp) for the
      /// moment..
      this->transport_jacobian() = D * this->transport_jacobian();
    }

    DETRAY_HOST_DEVICE
    constexpr vector3_type dtds() const { return {0.f, 0.f, 0.f}; }

    DETRAY_HOST_DEVICE
    constexpr scalar_type dqopds(
        const material<scalar_type>* /*unused*/) const {
      return 0.f;
    }
  };

  /// Take a step, regulared by a constrained step
  ///
  /// @param dist_to_next The straight line distance to the next surface
  /// @param stepping The state object of a stepper
  /// @param cfg The stepping configuration
  ///
  /// @returns returning the heartbeat, indicating if the stepping is alive
  DETRAY_HOST_DEVICE bool step(
      const scalar_type dist_to_next, state& stepping,
      const stepping::config& cfg, const bool = true,
      const material<scalar_type>* /*unused*/ = nullptr) const {
    // In case of an overlap do nothing
    if (math::fabs(dist_to_next) <= 1.f * unit<scalar_type>::um) {
      stepping.run_inspector(cfg, "Step skipped (Overlap): ", dist_to_next);
      return true;
    }

    // Straight line stepping: The distance given by the navigator is exact
    stepping.set_step_size(dist_to_next);

    DETRAY_VERBOSE_HOST_DEVICE("Take step: %f mm", stepping.step_size());

    // Check constraints
    if (const scalar_type max_step =
            stepping.constraints().template size<>(stepping.direction());
        math::fabs(stepping.step_size()) > math::fabs(max_step)) {
      // Run inspection before step size is cut
      stepping.run_inspector(cfg, "Before constraint: ", dist_to_next);

      stepping.set_step_size(max_step);
    }

    // Update track state
    stepping.advance_track();

    // Advance jacobian transport
    if (cfg.do_covariance_transport) {
      stepping.advance_jacobian();
    }

    // Count the number of steps
    stepping.count_trials();

    // Run inspection if needed
    stepping.run_inspector(cfg, "Step complete: ", dist_to_next);

    return true;
  }
};

}  // namespace detray
