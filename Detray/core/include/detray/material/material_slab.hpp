// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/material/material.hpp"
#include "detray/material/predefined_materials.hpp"

// System include(s)
#include <limits>
#include <ostream>

namespace detray {

// Slab structure to be mapped on the mask (plane, cylinder)
template <concepts::scalar scalar_t>
struct material_slab {
  using scalar_type = scalar_t;
  using material_type = material<scalar_t>;

  constexpr material_slab() = default;

  /// Constructor
  /// @param material is the elemental or mixture material
  /// @param thickness is the thickness of the slab
  constexpr material_slab(const material_type& material, scalar_type thickness)
      : m_material(material),
        m_thickness(thickness),
        m_thickness_in_X0(thickness / material.X0()),
        m_thickness_in_L0(thickness / material.L0()) {}

  /// Equality operator
  ///
  /// @param rhs is the right hand side to be compared to
  DETRAY_HOST_DEVICE
  constexpr bool operator==(const material_slab& rhs) const {
    return (m_material == rhs.get_material() && m_thickness == rhs.thickness());
  }

  /// Boolean operator
  DETRAY_HOST_DEVICE
  constexpr explicit operator bool() const {
    if (m_thickness <= std::numeric_limits<scalar_type>::epsilon() ||
        m_thickness == std::numeric_limits<scalar_type>::max() ||
        m_material == vacuum<scalar_type>() ||
        m_material.mass_density() == 0.f || m_material.molar_density() == 0.f) {
      return false;
    }
    return true;
  }

  /// Access the (average) material parameters.
  DETRAY_HOST_DEVICE
  constexpr const material_type& get_material() const { return m_material; }
  /// Return the thickness.
  DETRAY_HOST_DEVICE
  constexpr scalar_type thickness() const { return m_thickness; }
  /// Return the radiation length fraction.
  DETRAY_HOST_DEVICE
  constexpr scalar_type thickness_in_X0() const { return m_thickness_in_X0; }
  /// Return the nuclear interaction length fraction.
  DETRAY_HOST_DEVICE
  constexpr scalar_type thickness_in_L0() const { return m_thickness_in_L0; }

  /// @returns the path segment through the material
  ///
  /// @param cos_inc_angle cosine of the track incidence angle
  DETRAY_HOST_DEVICE constexpr scalar_type path_segment(
      const scalar_type cos_inc_angle, const scalar_type = 0.f) const {
    return m_thickness / cos_inc_angle;
  }

  /// @returns the path segment through the material in X0
  ///
  /// @param cos_inc_angle cosine of the track incidence angle
  DETRAY_HOST_DEVICE constexpr scalar_type path_segment_in_X0(
      const scalar_type cos_inc_angle, const scalar_type = 0.f) const {
    return m_thickness_in_X0 / cos_inc_angle;
  }

  /// @returns the path segment through the material in L0
  ///
  /// @param cos_inc_angle cosine of the track incidence angle
  DETRAY_HOST_DEVICE constexpr scalar_type path_segment_in_L0(
      const scalar_type cos_inc_angle, const scalar_type = 0.f) const {
    return m_thickness_in_L0 / cos_inc_angle;
  }

  /// @returns a string stream that prints the material details
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& os, const material_slab& mat) {
    os << "slab: ";
    os << mat.get_material();
    os << " | thickness = " << mat.thickness() << " mm";

    return os;
  }

 private:
  material_type m_material = {};
  scalar_type m_thickness = std::numeric_limits<scalar_type>::epsilon();
  scalar_type m_thickness_in_X0 = std::numeric_limits<scalar_type>::epsilon();
  scalar_type m_thickness_in_L0 = std::numeric_limits<scalar_type>::epsilon();
};

}  // namespace detray
