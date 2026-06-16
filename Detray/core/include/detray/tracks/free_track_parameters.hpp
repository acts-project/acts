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
#include "detray/definitions/track_parametrization.hpp"

// System include(s)
#include <ostream>

namespace detray {

template <concepts::algebra algebra_t>
struct free_parameters_vector {
  /// @name Type definitions for the struct
  /// @{
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using point3_type = dpoint3D<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;

  // Shorthand vector type related to free track parameters.
  using vector_type = free_vector<algebra_t>;

  /// @}

  /// Default constructor
  free_parameters_vector() = default;

  /// Construct from a 6-dim vector of parameters
  DETRAY_HOST_DEVICE
  explicit free_parameters_vector(const vector_type& vec) : m_vector(vec) {
    assert(!this->is_invalid());
  }

  /// Construct from single parameters
  ///
  /// @param pos the global position
  /// @param time the time
  /// @param mom the global track momentum 3-vector
  /// @param q the particle charge
  DETRAY_HOST_DEVICE
  free_parameters_vector(const point3_type& pos, const scalar_type time,
                         const vector3_type& mom, const scalar_type q) {
    getter::set_block(m_vector, pos, e_free_pos0, 0u);
    getter::element(m_vector, e_free_time, 0u) = time;

    scalar_type p = vector::norm(mom);
    vector3_type mom_norm = vector::normalize(mom);
    getter::set_block(m_vector, mom_norm, e_free_dir0, 0u);
    getter::element(m_vector, e_free_qoverp, 0u) = q / p;

    assert(!this->is_invalid());
  }

  /// @param rhs is the left hand side params for comparison
  DETRAY_HOST_DEVICE
  bool operator==(const free_parameters_vector& rhs) const {
    for (unsigned int i = 0u; i < e_free_size; i++) {
      if (math::fabs((*this)[i] - rhs[i]) >
          std::numeric_limits<scalar_type>::epsilon()) {
        return false;
      }
    }

    return true;
  }

  /// Convenience access to the track parameters - const
  DETRAY_HOST_DEVICE
  scalar_type operator[](std::size_t i) const {
    return getter::element(m_vector, static_cast<unsigned int>(i), 0u);
  }

  /// Convenience access to the track parameters - non-const
  DETRAY_HOST_DEVICE
  decltype(auto) operator[](std::size_t i) {
    return getter::element(m_vector, static_cast<unsigned int>(i), 0u);
  }

  /// @returns the global track position
  DETRAY_HOST_DEVICE
  point3_type pos() const {
    return {getter::element(m_vector, e_free_pos0, 0u),
            getter::element(m_vector, e_free_pos1, 0u),
            getter::element(m_vector, e_free_pos2, 0u)};
  }

  /// Set the global track position
  DETRAY_HOST_DEVICE
  void set_pos(const vector3_type& pos) {
    assert(math::isfinite(pos[0]));
    assert(math::isfinite(pos[1]));
    assert(math::isfinite(pos[2]));
    getter::set_block(m_vector, pos, e_free_pos0, 0u);
  }

  /// @returns the normalized, global track direction
  DETRAY_HOST_DEVICE
  vector3_type dir() const {
    return {getter::element(m_vector, e_free_dir0, 0u),
            getter::element(m_vector, e_free_dir1, 0u),
            getter::element(m_vector, e_free_dir2, 0u)};
  }

  /// Set the global track direction
  /// @note Must be normalized!
  DETRAY_HOST_DEVICE
  void set_dir(const vector3_type& dir) {
    assert(math::isfinite(dir[0]));
    assert(math::isfinite(dir[1]));
    assert(math::isfinite(dir[2]));
    assert(algebra::approx_equal(vector::norm(dir), scalar_type{1},
                                 scalar_type{1e-5f}));
    getter::set_block(m_vector, dir, e_free_dir0, 0u);
  }

  /// @returns the time
  DETRAY_HOST_DEVICE
  scalar_type time() const {
    return getter::element(m_vector, e_free_time, 0u);
  }

  /// Set the time
  DETRAY_HOST_DEVICE
  void set_time(const scalar_type t) {
    assert(math::isfinite(t));
    getter::element(m_vector, e_free_time, 0u) = t;
  }

  /// @returns the q/p value
  DETRAY_HOST_DEVICE
  scalar_type qop() const {
    return getter::element(m_vector, e_free_qoverp, 0u);
  }

  /// Set the q/p value
  DETRAY_HOST_DEVICE
  void set_qop(const scalar_type qop) {
    assert(math::isfinite(qop));
    getter::element(m_vector, e_free_qoverp, 0u) = qop;
  }

  /// @returns the q/p_T value
  DETRAY_HOST_DEVICE
  scalar_type qopT() const {
    const vector3_type dir = this->dir();
    assert(vector::perp(dir) != 0.f);
    return getter::element(m_vector, e_free_qoverp, 0u) / vector::perp(dir);
  }

  /// @returns the q/p_z value
  DETRAY_HOST_DEVICE
  scalar_type qopz() const {
    const vector3_type dir = this->dir();
    return getter::element(m_vector, e_free_qoverp, 0u) / dir[2];
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
    assert(this->qop() != 0.f);
    assert(q * qop() > 0.f);
    return math::fabs(q / this->qop() * vector::perp(this->dir()));
  }

  /// @returns the absolute momentum z-component
  DETRAY_HOST_DEVICE
  scalar_type pz(const scalar_type q) const {
    assert(this->qop() != 0.f);
    assert(q * qop() > 0.f);
    return math::fabs(q / this->qop() * this->dir()[2]);
  }

  /// @returns true if the parameter vector contains invalid elements or the
  /// direction is not normalized
  DETRAY_HOST_DEVICE
  constexpr bool is_invalid() const {
    bool inv_elem{false};
    for (std::size_t i = 0u; i < e_free_size; ++i) {
      inv_elem = inv_elem || !math::isfinite(getter::element(m_vector, i, 0u));
    }
    return inv_elem ||
           !algebra::approx_equal(vector::norm(dir()), scalar_type{1},
                                  scalar_type{1e-5f});
  }

 private:
  /// Transform to a string for debugging output
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& out_stream,
                                  const free_parameters_vector& fparam) {
    out_stream << fparam.m_vector;
    return out_stream;
  }

  vector_type m_vector = matrix::zero<vector_type>();
};

/// The free track parameters consist only of the parameter vector
template <concepts::algebra algebra_t>
using free_track_parameters = free_parameters_vector<algebra_t>;

}  // namespace detray
