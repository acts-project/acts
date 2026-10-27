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
#include "detray/material/material.hpp"
#include "detray/material/predefined_materials.hpp"

// System include(s)
#include <limits>
#include <ostream>

namespace detray {

/// Slab material attached to a planar or cylindrical surface mask.
template <concepts::scalar scalar_t>
struct material_slab {
  /// Scalar type used for thickness and material properties.
  using scalar_type = scalar_t;
  /// Bulk material type represented by this slab.
  using material_type = material<scalar_t>;

  constexpr material_slab() = default;

  /// Construct a slab from a material and its thickness.
  /// @param material bulk or mixture material in the slab
  /// @param thickness geometric thickness of the slab
  constexpr material_slab(const material_type& material, scalar_type thickness)
      : m_material(material),
        m_thickness(thickness),
        m_thickness_in_X0(thickness / material.X0()),
        m_thickness_in_L0(thickness / material.L0()) {}

  /// Compare the slab material and thickness.
  /// @param rhs slab to compare with this value
  /// @return whether material and thickness match
  DETRAY_HOST_DEVICE
  constexpr bool operator==(const material_slab& rhs) const {
    return (m_material == rhs.get_material() && m_thickness == rhs.thickness());
  }

  /// Test whether the slab represents non-vacuum material with positive
  /// thickness.
  /// @return true when the slab contributes material
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

  /// Return the material properties assigned to the slab.
  /// @return slab material
  DETRAY_HOST_DEVICE
  constexpr const material_type& get_material() const { return m_material; }
  /// Return the geometric slab thickness.
  /// @return thickness in length units
  DETRAY_HOST_DEVICE
  constexpr scalar_type thickness() const { return m_thickness; }
  /// Return thickness divided by the radiation length.
  /// @return fraction of one radiation length
  DETRAY_HOST_DEVICE
  constexpr scalar_type thickness_in_X0() const { return m_thickness_in_X0; }
  /// Return thickness divided by the nuclear interaction length.
  /// @return fraction of one nuclear interaction length
  DETRAY_HOST_DEVICE
  constexpr scalar_type thickness_in_L0() const { return m_thickness_in_L0; }

  /// Compute the path length through the slab at an incidence angle.
  /// @return path length through the slab
  ///
  /// @param cos_inc_angle cosine of the track incidence angle
  DETRAY_HOST_DEVICE constexpr scalar_type path_segment(
      const scalar_type cos_inc_angle,
      const scalar_type /*unused*/ = 0.f) const {
    return m_thickness / cos_inc_angle;
  }

  /// Compute the path length through the slab at an incidence angle.
  /// @return path length through the slab in X0
  ///
  /// @param cos_inc_angle cosine of the track incidence angle
  DETRAY_HOST_DEVICE constexpr scalar_type path_segment_in_X0(
      const scalar_type cos_inc_angle,
      const scalar_type /*unused*/ = 0.f) const {
    return m_thickness_in_X0 / cos_inc_angle;
  }

  /// Compute the path length through the slab at an incidence angle.
  /// @return path length through the slab in L0
  ///
  /// @param cos_inc_angle cosine of the track incidence angle
  DETRAY_HOST_DEVICE constexpr scalar_type path_segment_in_L0(
      const scalar_type cos_inc_angle,
      const scalar_type /*unused*/ = 0.f) const {
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
