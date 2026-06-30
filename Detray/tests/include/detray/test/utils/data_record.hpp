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
#include "detray/definitions/navigation.hpp"
#include "detray/propagator/actors/parameter_updater.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/tracks/free_track_parameters.hpp"

// System include(s)
#include <optional>

namespace detray {

/// Record of a surface intersection along a track
template <typename detector_t>
struct intersection_record {
  using algebra_type = typename detector_t::algebra_type;
  using scalar_type = dscalar<algebra_type>;
  using point3_type = dpoint3D<algebra_type>;
  using vector3_type = dvector3D<algebra_type>;
  using intersection_type =
      intersection2D<typename detector_t::surface_type, algebra_type,
                     intersection::contains_pos>;

  constexpr intersection_record() = default;

  /// The particle charge is not known in the navigation, but might be
  /// provided in a different context
  DETRAY_HOST_DEVICE
  constexpr intersection_record(
      const point3_type& position, const vector3_type& direction,
      const intersection_type& intr,
      const scalar_type q = detray::detail::invalid_value<scalar_type>(),
      const scalar_type p = detray::detail::invalid_value<scalar_type>(),
      dindex v_idx = dindex_invalid)
      : pos{position},
        dir{direction},
        intersection{intr},
        vol_idx{v_idx},
        charge{q},
        p_mag{p} {}

  /// Current global track position
  point3_type pos{0.f, 0.f, 0.f};
  /// Current global track direction
  vector3_type dir{0.f, 0.f, 1.f};
  /// The intersection result
  intersection_type intersection{};
  /// Index of the volume the intersection was found in
  dindex vol_idx{dindex_invalid};
  /// Charge hypothesis of the particle (invalid value if not known)
  scalar_type charge{detray::detail::invalid_value<scalar_type>()};
  /// Current momentum magnitude of the particle
  scalar_type p_mag{1.f};

  /// @returns the information as free track parameters
  DETRAY_HOST_DEVICE
  free_track_parameters<algebra_type> track_param() const {
    return {pos, detail::invalid_value<scalar_type>(), dir, p_mag / charge};
  }
};

/// Data for a single step
template <concepts::algebra algebra_t>
struct step_record {
  using scalar_type = dscalar<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using track_param_type = free_track_parameters<algebra_t>;
  using bound_param_type = bound_track_parameters<algebra_t>;
  using free_matrix_type = free_matrix<algebra_t>;

  scalar_type step_size{0.f};
  scalar_type path_length{0.f};
  std::size_t n_total_trials{0u};
  navigation::direction nav_dir = navigation::direction::e_forward;
  geometry::identifier identifier{};
  track_param_type free_params{};
  bound_param_type bound_params{};
  free_matrix_type jacobian{};
};

/// Contains the material parameters and the pathlength
template <concepts::scalar scalar_t>
struct material_record {
  /// The surface the material belongs to
  geometry::identifier geo_id{};
  /// Pathlength of the track through the material
  scalar_t path{detail::invalid_value<scalar_t>()};
  /// Material thickness/radius
  scalar_t thickness{detail::invalid_value<scalar_t>()};
  /// Radiation length
  scalar_t mat_X0{0.f};
  /// Interaction length
  scalar_t mat_L0{0.f};
};

/// Data for a single step
template <typename detector_t>
struct propagation_record {
  using algebra_type = typename detector_t::algebra_type;
  using scalar_type = dscalar<algebra_type>;
  using free_param_type = free_track_parameters<algebra_type>;
  using bound_param_type = bound_track_parameters<algebra_type>;
  using free_matrix_type = free_matrix<algebra_type>;
  using intersection_type = intersection_record<detector_t>::intersection_type;

  /// @returns the free track parameters
  DETRAY_HOST_DEVICE
  scalar_type step_size() const { return m_step_record.step_size; }
  /// Set the new step size @param step_size
  DETRAY_HOST_DEVICE
  void step_size(scalar_type step_size) { m_step_record.step_size = step_size; }

  /// @returns the free track parameters
  DETRAY_HOST_DEVICE
  scalar_type path_length() const { return m_step_record.path_length; }
  /// Set the new path length @param path_length
  DETRAY_HOST_DEVICE
  void path_length(scalar_type path_length) {
    m_step_record.path_length = path_length;
  }

  /// @returns the free track parameters
  DETRAY_HOST_DEVICE
  std::size_t n_rkn_trials() const { return m_step_record.n_total_trials; }
  /// Set the new number of RKN iterations @param n_rkn_trials
  DETRAY_HOST_DEVICE
  void n_rkn_trials(std::size_t n_rkn_trials) {
    m_step_record.n_total_trials = n_rkn_trials;
  }

  /// @returns the free track parameters
  DETRAY_HOST_DEVICE
  scalar_type charge() const { return m_intersection_record.charge; }
  /// Set the new charge @param charge
  DETRAY_HOST_DEVICE
  void charge(scalar_type charge) { m_intersection_record.charge = charge; }

  /// @returns the current total momentum magnitude
  DETRAY_HOST_DEVICE
  scalar_type p_mag() const { return m_intersection_record.p_mag; }
  /// Set the new total momentum magnitude @param p_mag
  DETRAY_HOST_DEVICE
  void p_mag(scalar_type p_mag) { m_intersection_record.p_mag = p_mag; }

  /// @returns the current geometry identifier if the step landed on surface
  DETRAY_HOST_DEVICE
  geometry::identifier identifier() const {
    assert(intersection().is_outside() &&
           m_step_record.bound_params.surface_link() ==
               intersection().surface().identifier());
    return m_step_record.bound_params.surface_link();
  }

  /// Set the new geometry identifier @param geo_id
  DETRAY_HOST_DEVICE
  void identifier(geometry::identifier geo_id) {
    m_step_record.bound_params.surface_link() = geo_id;
  }

  /// @returns the intersection data on the current surface (if any) - const
  DETRAY_HOST_DEVICE
  const intersection_type& intersection() const {
    return m_intersection_record.intersection;
  }
  /// @returns the intersection data on the current surface (if any)
  DETRAY_HOST_DEVICE
  intersection_type& intersection() {
    return m_intersection_record.intersection;
  }
  /// Set the new intersection data @param intr
  DETRAY_HOST_DEVICE
  void intersection(const intersection_type& intr) {
    m_intersection_record.intersection = intr;
  }

  /// @returns the free track parameters
  DETRAY_HOST_DEVICE
  free_param_type free_track_param() const { return m_step_record.free_params; }
  /// Set the new free track parameters @param free_param
  DETRAY_HOST_DEVICE
  void free_track_param(const free_param_type& free_param) {
    m_step_record.free_params = free_param;
  }

  /// @returns the bound track parameters at the current surface (if any) - const
  DETRAY_HOST_DEVICE
  bound_param_type& bound_track_param() { return m_step_record.bound_params; }
  /// Set the new bound parameters @param bound_param
  DETRAY_HOST_DEVICE
  void bound_track_param(const bound_param_type& bound_param) {
    m_step_record.bound_params = bound_param;
  }
  /// @returns the bound track parameters at the current surface (if any)
  DETRAY_HOST_DEVICE
  const bound_param_type& bound_track_param() const {
    return m_step_record.bound_params;
  }

  /// @returns transport Jacobian since the last surface - const
  DETRAY_HOST_DEVICE
  free_matrix_type& jacobian() { return m_step_record.jacobian; }
  /// @returns transport Jacobian since the last surface
  DETRAY_HOST_DEVICE
  const free_matrix_type& jacobian() const { return m_step_record.jacobian; }
  /// Set the new Jacobian @param jac
  DETRAY_HOST_DEVICE
  void jacobian(const free_matrix_type& jac) { m_step_record.jacobian = jac; }

 private:
  /// Data that is collected at every step
  step_record<algebra_type> m_step_record;
  /// Data that is collected when the propagation reached a surface
  std::optional<intersection_record<detector_t>> m_intersection_record{};
  /// Data that is collected when the propagation reached a material surface
  std::optional<material_record<scalar_type>> m_material_record{};
};

}  // namespace detray
