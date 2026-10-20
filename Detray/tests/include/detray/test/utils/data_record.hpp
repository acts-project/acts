// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/core/concepts.hpp"
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
template <concepts::detector detector_t>
struct intersection_record {
 private:
  using value_t = dvalue<typename detector_t::algebra_type>;

 public:
  using algebra_type = dalgebra<typename detector_t::algebra_type, double>;
  using scalar_type = dscalar<algebra_type>;
  using point3_type = dpoint3D<algebra_type>;
  using vector3_type = dvector3D<algebra_type>;
  using intersection_type = intersection2D<typename detector_t::surface_type,
                                           typename detector_t::algebra_type,
                                           intersection::contains_pos>;

  constexpr intersection_record() = default;

  /// The particle charge is not known in the navigation, but might be
  /// provided in a different context
  DETRAY_HOST_DEVICE
  constexpr intersection_record(
      const typename detector_t::point3_type& position,
      const typename detector_t::vector3_type& direction,
      const intersection_type& intr,
      const double q = detray::detail::invalid_value<double>(),
      const double p = detray::detail::invalid_value<double>(),
      const dindex v_idx = dindex_invalid)
      : m_pos{algebra::cast_to<double>(position)},
        m_dir{algebra::cast_to<double>(direction)},
        m_intersection{intr},
        m_vol_idx{v_idx},
        m_charge{q},
        m_p_mag{p} {}

  /// Getters
  /// @{
  template <concepts::value V = value_t>
  DETRAY_HOST_DEVICE auto position() const {
    return algebra::cast_to<V>(m_pos);
  }
  template <concepts::value V = value_t>
  DETRAY_HOST_DEVICE auto direction() const {
    return algebra::cast_to<V>(m_dir);
  }
  DETRAY_HOST_DEVICE
  const intersection_type& intersection() const { return m_intersection; }
  DETRAY_HOST_DEVICE
  intersection_type& intersection() { return m_intersection; }
  DETRAY_HOST_DEVICE
  dindex volume() const { return m_vol_idx; }
  template <concepts::value V = value_t>
  DETRAY_HOST_DEVICE auto charge() const {
    return algebra::cast_to<V>(m_charge);
  }
  template <concepts::value V = value_t>
  DETRAY_HOST_DEVICE auto p() const {
    return algebra::cast_to<V>(m_p_mag);
  }

  /// @returns the information as free track parameters
  DETRAY_HOST_DEVICE
  free_track_parameters<typename detector_t::algebra_type> track_param() const {
    assert(m_p_mag > 0.f);
    return {algebra::cast_to<value_t>(m_pos),
            detail::invalid_value<typename detector_t::scalar_type>(),
            algebra::cast_to<value_t>(m_p_mag * m_dir),
            algebra::cast_to<value_t>(m_charge)};
  }
  /// @}

  /// Setters
  /// @{
  template <concepts::point3D point_t>
  DETRAY_HOST_DEVICE void position(const point_t& position) {
    m_pos = algebra::cast_to<double>(position);
  }
  template <concepts::vector3D vector_t>
  DETRAY_HOST_DEVICE void direction(const vector_t& dir) {
    m_dir = algebra::cast_to<double>(dir);
  }
  DETRAY_HOST_DEVICE
  void intersection(const intersection_type& intr) { m_intersection = intr; }
  DETRAY_HOST_DEVICE
  void volume(dindex vol_idx) { m_vol_idx = vol_idx; }

  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE void charge(const scalar_t q) {
    m_charge = algebra::cast_to<double>(q);
  }
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE void p(const scalar_t p) {
    m_p_mag = algebra::cast_to<double>(p);
  }
  /// @}

 private:
  /// Current global track position
  point3_type m_pos{0., 0., 0.};
  /// Current global track direction
  vector3_type m_dir{0., 0., 1.};
  /// The intersection result
  intersection_type m_intersection{};
  /// Index of the volume the intersection was found in
  dindex m_vol_idx{dindex_invalid};
  /// Charge hypothesis of the particle (invalid value if not known)
  scalar_type m_charge{detray::detail::invalid_value<scalar_type>()};
  /// Current momentum magnitude of the particle
  scalar_type m_p_mag{1. * unit<scalar_type>::GeV};
};

/// Data for a single step
template <concepts::algebra algebra_t>
struct step_record {
  using algebra_type = dalgebra<algebra_t, double>;
  using value_t = dvalue<algebra_t>;
  using scalar_type = dscalar<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using track_param_type = free_track_parameters<algebra_t>;
  using bound_param_type = bound_track_parameters<algebra_t>;
  using free_matrix_type = free_matrix<algebra_t>;

  scalar_type m_step_size{0.};
  scalar_type m_path_length{0.};
  std::size_t m_n_total_trials{0u};
  navigation::direction m_nav_dir = navigation::direction::e_forward;
  geometry::identifier m_identifier{};
  track_param_type m_free_params{};
  bound_param_type m_bound_params{};
  free_matrix_type m_jacobian{};

  /// Getters
  /// @{
  template <concepts::value V = value_t>
  DETRAY_HOST_DEVICE scalar_type step_size() const {
    return algebra::cast_to<V>(m_step_size);
  }
  template <concepts::value V = value_t>
  DETRAY_HOST_DEVICE scalar_type path_length() const {
    return algebra::cast_to<V>(m_path_length);
  }
  DETRAY_HOST_DEVICE std::size_t n_trials() const { return m_n_total_trials; }
  DETRAY_HOST_DEVICE navigation::direction nav_direction() const {
    return m_nav_dir;
  }
  DETRAY_HOST_DEVICE geometry::identifier identifier() const {
    return m_identifier;
  }
  template <concepts::value V = value_t>
  DETRAY_HOST_DEVICE track_param_type free_param() const {
    return algebra::cast_to<V>(m_free_params);
  }
  template <concepts::value V = value_t>
  DETRAY_HOST_DEVICE bound_param_type bound_param() const {
    return algebra::cast_to<V>(m_bound_params);
  }
  template <concepts::value V = value_t>
  DETRAY_HOST_DEVICE free_matrix_type jacobian() const {
    return algebra::cast_to<V>(m_jacobian);
  }
  /// @}

  /// Setters
  /// @{
  DETRAY_HOST_DEVICE void step_size(const scalar_type step_size) {
    m_step_size = algebra::cast_to<double>(step_size);
  }
  DETRAY_HOST_DEVICE void path_length(const scalar_type path_length) {
    m_path_length = algebra::cast_to<double>(path_length);
  }
  DETRAY_HOST_DEVICE void n_trials(const std::size_t n_trial) {
    m_n_total_trials = n_trial;
  }
  DETRAY_HOST_DEVICE void nav_direction(const navigation::direction nav_dir) {
    m_nav_dir = nav_dir;
  }
  DETRAY_HOST_DEVICE void identifier(const geometry::identifier geo_id) {
    m_identifier = geo_id;
  }
  template <typename free_param_t>
  DETRAY_HOST_DEVICE void free_param(const free_param_t& param) {
    m_free_params = algebra::cast_to<double>(param);
  }
  template <typename bound_param_t>
  DETRAY_HOST_DEVICE void bound_param(const bound_param_t& param) {
    m_bound_params = algebra::cast_to<double>(param);
  }
  template <typename jacobian_t>
  DETRAY_HOST_DEVICE void jacobian(const jacobian_t& jac) {
    m_jacobian = algebra::cast_to<double>(jac);
  }
  /// @}
};

/// Contains the material parameters and the pathlength
struct material_record {
  /// The surface the material belongs to
  geometry::identifier geo_id{};
  /// Pathlength of the track through the material
  double path{detail::invalid_value<double>()};
  /// Material thickness/radius
  double thickness{detail::invalid_value<double>()};
  /// Radiation length
  double mat_X0{0.};
  /// Interaction length
  double mat_L0{0.};
};

/// Data for a single step
template <concepts::detector detector_t>
struct propagation_record {
  using algebra_type = typename detector_t::algebra_type;
  using value_type = dvalue<algebra_type>;
  using scalar_type = dscalar<algebra_type>;
  using step_record_type = step_record<algebra_type>;
  using free_param_type = typename step_record_type::track_param_type;
  using bound_param_type = typename step_record_type::bound_param_type;
  using free_matrix_type = typename step_record_type::free_matrix_type;
  using intersection_type = intersection_record<detector_t>::intersection_type;

  /// @returns the free track parameters
  DETRAY_HOST_DEVICE
  scalar_type step_size() const { return m_step_record.step_size(); }
  /// Set the new step size @param step_size
  DETRAY_HOST_DEVICE
  void step_size(scalar_type step_size) { m_step_record.step_size(step_size); }

  /// @returns the free track parameters
  DETRAY_HOST_DEVICE
  scalar_type path_length() const { return m_step_record.path_length(); }
  /// Set the new path length @param path_length
  DETRAY_HOST_DEVICE
  void path_length(scalar_type path_length) {
    m_step_record.path_length(path_length);
  }

  /// @returns the free track parameters
  DETRAY_HOST_DEVICE
  std::size_t n_rkn_trials() const { return m_step_record.n_trials(); }
  /// Set the new number of RKN iterations @param n_rkn_trials
  DETRAY_HOST_DEVICE
  void n_rkn_trials(std::size_t n_rkn_trials) {
    m_step_record.n_trials(n_rkn_trials);
  }

  /// @returns the free track parameters
  DETRAY_HOST_DEVICE
  scalar_type charge() const {
    return m_intersection_record.template charge<value_type>();
  }
  /// Set the new charge @param charge
  DETRAY_HOST_DEVICE
  void charge(scalar_type charge) { m_intersection_record.charge(charge); }

  /// @returns the current total momentum magnitude
  DETRAY_HOST_DEVICE
  scalar_type p_mag() const {
    return m_intersection_record.template p<value_type>();
  }
  /// Set the new total momentum magnitude @param p_mag
  DETRAY_HOST_DEVICE
  void p_mag(scalar_type p_mag) { m_intersection_record.p(p_mag); }

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
    m_step_record.bound_param().surface_link() = geo_id;
  }

  /// @returns the intersection data on the current surface (if any) - const
  DETRAY_HOST_DEVICE
  const intersection_type& intersection() const {
    return m_intersection_record.intersection();
  }
  /// @returns the intersection data on the current surface (if any)
  DETRAY_HOST_DEVICE
  intersection_type& intersection() {
    return m_intersection_record.intersection();
  }
  /// Set the new intersection data @param intr
  DETRAY_HOST_DEVICE
  void intersection(const intersection_type& intr) {
    m_intersection_record.intersection(intr);
  }

  /// @returns the free track parameters
  DETRAY_HOST_DEVICE
  free_param_type free_track_param() const {
    return m_step_record.free_param();
  }
  /// Set the new free track parameters @param free_param
  DETRAY_HOST_DEVICE
  void free_track_param(const free_param_type& free_param) {
    m_step_record.free_param(free_param);
  }

  /// @returns the bound track parameters at the current surface (if any) - const
  DETRAY_HOST_DEVICE
  bound_param_type& bound_track_param() { return m_step_record.bound_param(); }
  /// Set the new bound parameters @param bound_param
  DETRAY_HOST_DEVICE
  void bound_track_param(const bound_param_type& bound_param) {
    m_step_record.bound_param(bound_param);
  }
  /// @returns the bound track parameters at the current surface (if any)
  DETRAY_HOST_DEVICE
  const bound_param_type& bound_track_param() const {
    return m_step_record.bound_params;
  }

  /// @returns transport Jacobian since the last surface - const
  DETRAY_HOST_DEVICE
  free_matrix_type& jacobian() { return m_step_record.jacobian(); }
  /// @returns transport Jacobian since the last surface
  DETRAY_HOST_DEVICE
  const free_matrix_type& jacobian() const { return m_step_record.jacobian(); }
  /// Set the new Jacobian @param jac
  DETRAY_HOST_DEVICE
  void jacobian(const free_matrix_type& jac) { m_step_record.jacobian(jac); }

 private:
  /// Data that is collected at every step
  step_record<algebra_type> m_step_record;
  /// Data that is collected when the propagation reached a surface
  std::optional<intersection_record<detector_t>> m_intersection_record{};
  /// Data that is collected when the propagation reached a material surface
  std::optional<material_record> m_material_record{};
};

}  // namespace detray
