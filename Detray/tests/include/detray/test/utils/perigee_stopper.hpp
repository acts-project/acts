// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/line.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/constrained_step.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/curvilinear_frame.hpp"

namespace detray {

/// Actor that exits from a navigation stream, if it has reached the perigee
/// @note Currently only works with step constraints (will be changed in the
/// future)
template <concepts::algebra algebra_t>
struct perigee_stopper : public base_actor {
  using scalar_t = dscalar<algebra_t>;

  struct state {
    /// Radius around the beamline, where to test against the perigee
    /// Outside this radius, the track is considered as not originating
    /// from the IP
    scalar_t m_stopping_radius{10.f * unit<scalar_t>::mm};
    /// Tolerance under which to consider the track at the perigee
    // @TODO Make smaller once overstepping is solved
    scalar_t m_on_perigee_tol{100.f * unit<scalar_t>::um};
    /// Index of the innermost volume for this detector: Convention is 0
    unsigned int m_inner_vol_idx{0u};
  };

  /// Intersects a linear track approximation with the perigee and exits
  /// the navigation, if the perigee is reached.
  ///
  /// @param prop_state state of the propagation
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE void operator()(state &actor_state,
                                     propagator_state_t &prop_state) const {
    using detector_t = typename propagator_state_t::detector_type;
    using perigee_intersector_t = ray_intersector<line_circular, algebra_t>;

    // Nothing left to do. Propagation will exit successfully on its own
    auto &navigation = prop_state.navigation();
    if (navigation.finished()) {
      return;
    }

    // Only check this in the innermost volume and during backward nav.
    if (navigation.volume() != actor_state.m_inner_vol_idx ||
        navigation.direction() != navigation::direction::e_backward) {
      return;
    }

    // Volume that contains the IP
    const tracking_volume inner_vol{navigation.detector(),
                                    actor_state.m_inner_vol_idx};

    // Stop at the perigee (cylindrical detectors only)
    if (inner_vol.id() != volume_id::e_cylinder) {
      return;
    }

    // At least the exit portal should be reachable
    if (navigation.cache_exhausted()) {
      navigation.abort("Perigee stopper has no next candidate");
      prop_state.heartbeat(false);
      return;
    }

    auto &stepping = prop_state.stepping();
    auto &track = stepping();

    // Linear track approximation
    const detail::ray<algebra_t> trk_approx{track.pos(), -1.f * track.dir()};

    // Check the stopping radius (outside the track will not be stopped)
    assert(actor_state.m_stopping_radius > 0.f);
    constexpr scalar_t max_hz{detail::invalid_value<scalar_t>()};

    using line_mask_t = mask<line_circular, algebra_t>;
    const line_mask_t perigee_mask{
        static_cast<typename line_mask_t::links_type>(
            actor_state.m_inner_vol_idx),
        actor_state.m_stopping_radius, max_hz};

    // The perigee is not linked to any detector surface
    constexpr typename detector_t::surface_type inv_sf{};
    const dtransform3D<algebra_t> identity{};
    constexpr scalar_t mask_tolerance{0.f};
    constexpr scalar_t overstep_tolerance{-detail::invalid_value<scalar_t>()};

    const auto perigee_intr =
        perigee_intersector_t{}(trk_approx, inv_sf, perigee_mask, identity,
                                mask_tolerance, overstep_tolerance);

    scalar_t dist_to_cand{std::as_const(navigation).target().path()};
    if (perigee_intr.is_probably_inside() &&
        math::fabs(perigee_intr.path()) < math::fabs(dist_to_cand)) {
      assert(actor_state.m_on_perigee_tol > 0.f);

      // The track has reached the perigee: "exit success"
      if (math::fabs(perigee_intr.path()) <= actor_state.m_on_perigee_tol) {
        navigation.exit();
        prop_state.heartbeat(false);
      } else {
        // @TODO: Use a guided navigator for this in order to catch
        // overstepping correctly
        stepping.template set_constraint<step::constraint::e_actor>(
            perigee_intr.path());
      }
    }
  }
};

}  // namespace detray
