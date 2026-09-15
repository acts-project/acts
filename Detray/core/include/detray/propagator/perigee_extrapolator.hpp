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

#pragma once

// Project include(s).
#include "detray/definitions/detail/macros.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/shapes/line.hpp"
#include "detray/material/material_rod.hpp"
#include "detray/navigation/direct_navigator.hpp"
#include "detray/navigation/external_surface.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/concepts.hpp"
#include "detray/propagator/propagation_config.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/curvilinear_frame.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"
#include "detray/utils/ranges/single.hpp"

namespace detray {

/// Templated propagator class, using a stepper and a navigator object in
/// succession.
///
/// @tparam stepper_t for the transport
/// @tparam navigator_t for the navigation
template <typename detector_t, typename stepper_t,
          concepts::actor_chain actor_chain_t = actor_chain<>>
class perigee_extrapolator {
  using nav_link_t = typename detector_t::surface_type::navigation_link;

 public:
  using detector_type = detector_t;
  using algebra_type = typename detector_type::algebra_type;
  using scalar_type = dscalar<algebra_type>;

  using perigee_surface =
      external_surface<algebra_type, line_circular, material_rod<scalar_type>,
                       nav_link_t>;
  using navigator_type =
      direct_navigator<detector_t,
                       detray::ranges::single_view<perigee_surface>>;
  using intersection_type = typename navigator_type::intersection_type;

  using stepper_type = stepper_t;
  using free_track_parameters_type =
      typename stepper_t::free_track_parameters_type;
  using bound_track_parameters_type =
      typename stepper_t::bound_track_parameters_type;

  using actor_chain_type = actor_chain_t;

  using propagator_type =
      propagator<stepper_type, navigator_type, actor_chain_type>;

  using state = typename propagator_type::state;

  /// Construct from a propagator configuration
  DETRAY_HOST_DEVICE
  explicit constexpr perigee_extrapolator(const propagation::config &cfg)
      : m_propagator{cfg} {}

  /// Create the propagation state for the extrapolator
  ///
  /// @param track the free track parameters (initial state)
  /// @param det the detector
  ///
  /// @returns the propagation state.
  template <typename track_t>
  DETRAY_HOST_DEVICE state create_state(const track_t &track,
                                        const detector_t &det) const {
    return state{
        track, det,
        detray::ranges::single_view<perigee_surface>{m_perigee_surface},
        m_propagator.m_cfg.context};
  }

  /// Return whether or not the propagation completed successfully.
  ///
  /// @param propagation the state of the extrapolator
  ///
  /// @returns propagation success.
  DETRAY_HOST_DEVICE bool finished(const state &propagation) const {
    return m_propagator.finished(propagation);
  }

  /// Propagate method: Coordinates the calls of the stepper, navigator and
  /// all registered actors.
  ///
  /// @param propagation the state of a propagation flow
  /// @param actor_state_refs tule containing references to the actor states
  ///
  /// @returns propagation success.
  DETRAY_HOST_DEVICE auto extrapolate(
      state &propagation,
      typename actor_chain_type::state_ref_tuple actor_state_refs =
          typename actor_chain_type::state_ref_tuple{}) const {
    propagation.navigation().set_direction(navigation::direction::e_backward);
    m_propagator.propagate(propagation, actor_state_refs);

    assert(finished(propagation));
    assert(!propagation.stepping()().is_invalid());

    return detail::free_to_bound_vector<line2D<algebra_type>>(
        m_perigee_surface.transform(), propagation.stepping()());
  }

  /// Overload for empty actor chain
  DETRAY_HOST_DEVICE auto extrapolate(state &propagation) {
    // Run propagation
    propagation.navigation().set_direction(navigation::direction::e_backward);
    m_propagator.propagate(propagation);

    assert(finished(propagation));
    assert(!propagation.stepping()().is_invalid());

    return detail::free_to_bound_vector<line2D<algebra_type>>(
        m_perigee_surface.transform(), propagation.stepping()());
  }

 private:
  /// Material of the perigee surface
  static constexpr typename perigee_surface::material_type perigee_material{
      detray::vacuum<scalar_type>{}, std::numeric_limits<scalar_type>::max()};
  /// Line mask of the perigee
  static constexpr typename perigee_surface::mask_type perigee_mask{
      0u, std::numeric_limits<scalar_type>::max(),
      -std::numeric_limits<scalar_type>::max()};
  /// The perigee surface
  perigee_surface m_perigee_surface{dtransform3D<algebra_type>{}, perigee_mask,
                                    perigee_material, 0u,
                                    surface_id::e_passive};
  /// The propagator
  propagator_type m_propagator{};
};

}  // namespace detray
