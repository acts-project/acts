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
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/identifier.hpp"

// System include(s)
#include <ostream>

namespace detray {

template <concepts::algebra algebra_t>
struct bound_parameters_vector {
  /// @name Type definitions for the struct
  /// @{
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using point2_type = dpoint2D<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;

  // Underlying vector type related to bound track vector.
  using vector_type = bound_vector<algebra_t>;

  /// @}

  /// Default constructor
  bound_parameters_vector() = default;

  /// Construct from a 6-dim vector of parameters
  DETRAY_HOST_DEVICE
  explicit bound_parameters_vector(const vector_type& vec) : m_vector(vec) {
    assert(!this->is_invalid());
  }

  /// Construct from single parameters
  ///
  /// @param loc_p the bound local position
  /// @param phi the global phi angle of the track direction
  /// @param theta the global theta angle of the track direction
  /// @param qop the q/p value
  /// @param t the time
  DETRAY_HOST_DEVICE
  bound_parameters_vector(const point2_type& loc_p, const scalar_type phi,
                          const scalar_type theta, const scalar_type qop,
                          const scalar_type t) {
    getter::set_block(m_vector, loc_p, e_bound_loc0, 0u);
    getter::element(m_vector, e_bound_phi, 0u) = phi;
    getter::element(m_vector, e_bound_theta, 0u) = theta;
    getter::element(m_vector, e_bound_qoverp, 0u) = qop;
    getter::element(m_vector, e_bound_time, 0u) = t;

    assert(!this->is_invalid());
  }

  /// @param rhs is the left hand side params for comparison
  DETRAY_HOST_DEVICE
  bool operator==(const bound_parameters_vector& rhs) const {
    for (unsigned int i = 0u; i < e_bound_size; i++) {
      const auto lhs_val = getter::element(m_vector, i, 0u);
      const auto rhs_val = getter::element(rhs.vector(), i, 0u);

      if (math::fabs(lhs_val - rhs_val) >
          std::numeric_limits<scalar_type>::epsilon()) {
        return false;
      }
    }

    return true;
  }

  /// Convenience access to the track parameters - const
  DETRAY_HOST_DEVICE
  scalar_type operator[](const std::size_t i) const {
    return getter::element(m_vector, static_cast<unsigned int>(i), 0u);
  }

  /// Convenience access to the track parameters - non-const
  DETRAY_HOST_DEVICE
  decltype(auto) operator[](const std::size_t i) {
    return getter::element(m_vector, static_cast<unsigned int>(i), 0u);
  }

  /// Access the track parameters as a 6-dim vector - const
  DETRAY_HOST_DEVICE
  const vector_type& vector() const { return m_vector; }

  /// Access the track parameters as a 6-dim vector - non-const
  DETRAY_HOST_DEVICE
  vector_type& vector() { return m_vector; }

  /// Set the underlying vector
  DETRAY_HOST_DEVICE
  void set_vector(const vector_type& v,
                  [[maybe_unused]] const bool skip_check = false) {
    m_vector = v;
    assert(skip_check || !this->is_invalid());
  }

  /// @returns the bound local position
  DETRAY_HOST_DEVICE
  point2_type bound_local() const {
    return {getter::element(m_vector, e_bound_loc0, 0u),
            getter::element(m_vector, e_bound_loc1, 0u)};
  }

  /// Set the bound local position
  DETRAY_HOST_DEVICE
  void set_bound_local(const point2_type& pos) {
    assert(math::isfinite(pos[0]));
    assert(math::isfinite(pos[1]));
    getter::set_block(m_vector, pos, e_bound_loc0, 0u);
  }

  /// @returns the global phi angle
  DETRAY_HOST_DEVICE
  scalar_type phi() const { return getter::element(m_vector, e_bound_phi, 0u); }

  /// Set the global phi angle
  DETRAY_HOST_DEVICE
  void set_phi(const scalar_type phi) {
    assert(math::fabs(phi) <= constant<scalar_type>::pi);
    getter::element(m_vector, e_bound_phi, 0u) = phi;
  }

  /// @returns the global theta angle
  DETRAY_HOST_DEVICE
  scalar_type theta() const {
    return getter::element(m_vector, e_bound_theta, 0u);
  }

  /// Set the global theta angle
  DETRAY_HOST_DEVICE
  void set_theta(const scalar_type theta) {
    assert(0.f < theta);
    assert(theta <= constant<scalar_type>::pi);
    getter::element(m_vector, e_bound_theta, 0u) = theta;
  }

  /// @returns the global track direction
  DETRAY_HOST_DEVICE
  vector3_type dir() const {
    const scalar_type phi{getter::element(m_vector, e_bound_phi, 0u)};
    const scalar_type theta{getter::element(m_vector, e_bound_theta, 0u)};
    const scalar_type sinTheta{math::sin(theta)};

    return {math::cos(phi) * sinTheta, math::sin(phi) * sinTheta,
            math::cos(theta)};
  }

  /// @returns the time
  DETRAY_HOST_DEVICE
  scalar_type time() const {
    return getter::element(m_vector, e_bound_time, 0u);
  }

  /// Set the time
  DETRAY_HOST_DEVICE
  void set_time(const scalar_type t) {
    assert(math::isfinite(t));
    getter::element(m_vector, e_bound_time, 0u) = t;
  }

  /// @returns the q/p value
  DETRAY_HOST_DEVICE
  scalar_type qop() const {
    return getter::element(m_vector, e_bound_qoverp, 0u);
  }

  /// Set the q/p value
  DETRAY_HOST_DEVICE
  void set_qop(const scalar_type qop) {
    assert(math::isfinite(qop));
    getter::element(m_vector, e_bound_qoverp, 0u) = qop;
  }

  /// @returns the q/p_T value
  DETRAY_HOST_DEVICE
  scalar_type qopT() const {
    const scalar_type theta{getter::element(m_vector, e_bound_theta, 0u)};
    const scalar_type sinTheta{math::sin(theta)};
    assert(sinTheta != 0.f);
    return getter::element(m_vector, e_bound_qoverp, 0u) / sinTheta;
  }

  /// @returns the q/p_z value
  DETRAY_HOST_DEVICE
  scalar_type qopz() const {
    const scalar_type theta{getter::element(m_vector, e_bound_theta, 0u)};
    const scalar_type cosTheta{math::cos(theta)};
    assert(cosTheta != 0.f);
    return getter::element(m_vector, e_bound_qoverp, 0u) / cosTheta;
  }

  /// @returns the absolute momentum
  DETRAY_HOST_DEVICE
  scalar_type p(const scalar_type q) const {
    assert(qop() != 0.f);
    assert(q * qop() > 0.f);
    return q / qop();
  }

  /// @returns the global momentum 3-vector
  DETRAY_HOST_DEVICE
  vector3_type mom(const scalar_type q) const { return p(q) * dir(); }

  /// @returns the transverse momentum
  DETRAY_HOST_DEVICE
  scalar_type pT(const scalar_type q) const {
    assert(qop() != 0.f);
    assert(q * qop() > 0.f);
    return math::fabs(q / qop() * vector::perp(dir()));
  }

  /// @returns the absolute momentum z-component
  DETRAY_HOST_DEVICE
  scalar_type pz(const scalar_type q) const {
    assert(qop() != 0.f);
    assert(q * qop() > 0.f);
    return math::fabs(q / qop() * dir()[2]);
  }

  /// @param do_check toggle checking (e.g. don't trigger assertions for
  /// documented errors)
  /// @returns true if the parameter vector contains invalid elements
  DETRAY_HOST_DEVICE
  constexpr bool is_invalid(const bool do_check = true) const {
    if (!do_check) {
      return false;
    }
    bool inv_elem{false};
    bool is_all_zero{true};
    for (std::size_t i = 0u; i < e_bound_size; ++i) {
      inv_elem = inv_elem || !math::isfinite(getter::element(m_vector, i, 0u));
      is_all_zero =
          is_all_zero && (math::fabs(getter::element(m_vector, i, 0u)) == 0.f);
    }

    return (inv_elem || is_all_zero);
  }

 private:
  /// Transform to a string for debugging output
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& out_stream,
                                  const bound_parameters_vector& bparam) {
    out_stream << bparam.m_vector;
    return out_stream;
  }

  vector_type m_vector = matrix::zero<vector_type>();
};

/// Combine the bound track parameter vector with the covariance and associated
/// surface
template <concepts::algebra algebra_t>
struct bound_track_parameters : public bound_parameters_vector<algebra_t> {
  using base_type = bound_parameters_vector<algebra_t>;

  /// @name Type definitions for the struct
  /// @{
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using point2_type = dpoint2D<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;

  // Shorthand vector/matrix types related to bound track parameters.
  using parameter_vector_type = bound_parameters_vector<algebra_t>;
  using covariance_type = bound_matrix<algebra_t>;

  /// @}

  /// Default constructor sets the covaraicne to zero
  bound_track_parameters() = default;

  DETRAY_HOST_DEVICE
  bound_track_parameters(const geometry::identifier sf_idx,
                         const parameter_vector_type& vec,
                         const covariance_type& cov)
      : base_type(vec), m_covariance(cov), m_identifier(sf_idx) {}

  /// @param rhs is the left hand side params for comparison
  DETRAY_HOST_DEVICE
  bool operator==(const bound_track_parameters& rhs) const {
    if (m_identifier != rhs.surface_link()) {
      return false;
    }

    if (!base_type::operator==(rhs)) {
      return false;
    }

    for (unsigned int i = 0u; i < e_bound_size; i++) {
      for (unsigned int j = 0u; j < e_bound_size; j++) {
        const auto lhs_val = getter::element(m_covariance, i, j);
        const auto rhs_val = getter::element(rhs.covariance(), i, j);

        if (math::fabs(lhs_val - rhs_val) >
            std::numeric_limits<scalar_type>::epsilon()) {
          return false;
        }
      }
    }

    return true;
  }

  /// @returns the identifier of the associated surface
  DETRAY_HOST_DEVICE
  const geometry::identifier& surface_link() const { return m_identifier; }

  /// Set the identifier of the associated surface
  DETRAY_HOST_DEVICE
  void set_surface_link(geometry::identifier link) { m_identifier = link; }

  /// Set the track parameter vector
  DETRAY_HOST_DEVICE
  void set_parameter_vector(const parameter_vector_type& v,
                            const bool skip_check = false) {
    this->set_vector(v.vector(), skip_check);
  }

  /// @returns the track parameter covariance - non-const
  DETRAY_HOST_DEVICE
  covariance_type& covariance() { return m_covariance; }

  /// @returns the track parameter covariance - const
  DETRAY_HOST_DEVICE
  const covariance_type& covariance() const { return m_covariance; }

  /// Set the track parameter covariance
  DETRAY_HOST_DEVICE
  void set_covariance(const covariance_type& c) { m_covariance = c; }

  /// @param do_check toggle checking (e.g. don't trigger assertions for
  /// documented errors)
  /// @returns true if the parameter vector contains invalid elements
  DETRAY_HOST_DEVICE
  constexpr bool is_invalid(const bool do_check = true) const {
    if (!do_check) {
      return false;
    }
    if (base_type::is_invalid()) {
      return true;
    }

    // @TODO: Add tests positive semi-definite, check the determinant etc
    return (m_covariance == matrix::zero<covariance_type>());
  }

 private:
  /// Transform to a string for debugging output
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& out_stream,
                                  const bound_track_parameters& bparam) {
    out_stream << "Surface: " << bparam.m_identifier << std::endl;
    out_stream << "Param.:\n " << static_cast<parameter_vector_type>(bparam)
               << std::endl;
    out_stream << "Cov.:\n" << bparam.m_covariance;

    return out_stream;
  }

  /// Bound covaraicne matrix of the track parameters
  covariance_type m_covariance = matrix::zero<covariance_type>();
  /// Identifier of the surface the track parameters are associated with
  geometry::identifier m_identifier{};
};

}  // namespace detray
