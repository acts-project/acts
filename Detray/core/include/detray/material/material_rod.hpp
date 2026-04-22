// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/material/material.hpp"
#include "detray/material/predefined_materials.hpp"

// System include(s)
#include <limits>
#include <ostream>

namespace detray {

// Rod structure to be mapped on the line mask
template <concepts::scalar scalar_t>
struct material_rod {
  using scalar_type = scalar_t;
  using material_type = material<scalar_t>;

  constexpr material_rod() = default;

  constexpr material_rod(const material_type& material, scalar_type radius)
      : m_material(material),
        m_radius(radius),
        m_radius_in_X0(radius / material.X0()),
        m_radius_in_L0(radius / material.L0()) {}

  /// Equality operator
  ///
  /// @param rhs is the right hand side to be compared to
  DETRAY_HOST_DEVICE
  constexpr bool operator==(const material_rod<scalar_type>& rhs) const {
    return (m_material == rhs.get_material() && m_radius == rhs.radius());
  }

  /// Boolean operator
  DETRAY_HOST_DEVICE
  constexpr explicit operator bool() const {
    if (m_radius <= std::numeric_limits<scalar_type>::epsilon() ||
        m_material == vacuum<scalar_type>() || m_material.Z() == 0 ||
        m_material.mass_density() == 0) {
      return false;
    }
    return true;
  }

  /// Access the (average) material parameters.
  DETRAY_HOST_DEVICE
  constexpr const material_type& get_material() const { return m_material; }

  /// @returns the radius
  DETRAY_HOST_DEVICE
  constexpr scalar_type radius() const { return m_radius; }

  /// @returns the radius (alias function to keep interfaces homogeneous)
  DETRAY_HOST_DEVICE
  constexpr scalar_type thickness() const { return m_radius; }

  /// @returns the path segment through the material
  ///
  /// @param cos_inc_angle cosine of the track incidence angle
  /// @param approach the closes approach of the track to the wire
  DETRAY_HOST_DEVICE constexpr scalar_type path_segment(
      const scalar_type cos_inc_angle, const scalar_type approach) const {
    // Assume that is.local[0] is radial distance of line intersector
    if (math::fabs(approach) > m_radius) {
      return 0.f;
    }

    const scalar_type sin_inc_angle_2{1.f - cos_inc_angle * cos_inc_angle};

    return 2.f * math::sqrt((m_radius * m_radius - approach * approach) /
                            sin_inc_angle_2);
  }

  /// @returns the path segment through the material in X0
  ///
  /// @param cos_inc_angle cosine of the track incidence angle
  /// @param approach the closes approach of the track to the wire
  DETRAY_HOST_DEVICE constexpr scalar_type path_segment_in_X0(
      const scalar_type cos_inc_angle, const scalar_type approach) const {
    return this->path_segment(cos_inc_angle, approach) / m_material.X0();
  }

  /// @returns the path segment through the material in L0
  ///
  /// @param cos_inc_angle cosine of the track incidence angle
  /// @param approach the closes approach of the track to the wire
  DETRAY_HOST_DEVICE constexpr scalar_type path_segment_in_L0(
      const scalar_type cos_inc_angle, const scalar_type approach) const {
    return this->path_segment(cos_inc_angle, approach) / m_material.L0();
  }

  /// @returns a string stream that prints the material details
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& os, const material_rod& mat) {
    os << "rod: ";
    os << mat.get_material();
    os << " | radius = " << mat.radius() << "mm";

    return os;
  }

 private:
  material_type m_material = {};
  scalar_type m_radius = std::numeric_limits<scalar_type>::epsilon();
  scalar_type m_radius_in_X0 = std::numeric_limits<scalar_type>::epsilon();
  scalar_type m_radius_in_L0 = std::numeric_limits<scalar_type>::epsilon();
};

}  // namespace detray
