// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <cmath>
#include <limits>

#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @brief Internal representation of a space point for Acts seeding
/// It wraps an external space point representation and allows for shift
/// and variance setting
template <typename external_spacepoint_t>
class InternalSpacePoint {
 public:
  InternalSpacePoint() = delete;

  /// Constructor from arguments
  /// @param sp The external spacepoint object
  /// @param globalPos The represented global position
  /// @param offsetXY The beam/vertex offset in xy
  /// @param variance The space point variance
  InternalSpacePoint(const external_spacepoint_t& sp, const Vector3D& globalPos,
                     const Vector2D& offsetXY, const Vector2D& variance);

  InternalSpacePoint(const InternalSpacePoint<external_spacepoint_t>& sp);
  ~InternalSpacePoint() = default;

  InternalSpacePoint<external_spacepoint_t>& operator=(
      const InternalSpacePoint<external_spacepoint_t>&);

  const float& x() const { return m_x; }
  const float& y() const { return m_y; }
  const float& z() const { return m_z; }
  const float& radius() const { return m_r; }
  float phi() const { return atan2f(m_y, m_x); }
  const float& varianceR() const { return m_varianceR; }
  const float& varianceZ() const { return m_varianceZ; }
  const external_spacepoint_t& sp() const { return m_sp; }

 protected:
  float m_x;          /// x-coordinate in beam system coordinates
  float m_y;          /// y-coordinate in beam system coordinates
  float m_z;          /// z-coordinate in beam system coordinetes
  float m_r;          /// radius       in beam system coordinates
  float m_varianceR;  ///
  float m_varianceZ;  ///
  const external_spacepoint_t& m_sp;  /// external space point (referenced)
};

template <typename external_spacepoint_t>
inline InternalSpacePoint<external_spacepoint_t>::InternalSpacePoint(
    const external_spacepoint_t& sp, const Vector3D& globalPos,
    const Vector2D& offsetXY, const Vector2D& variance)
    : m_sp(sp) {
  m_x = globalPos.x() - offsetXY.x();
  m_y = globalPos.y() - offsetXY.y();
  m_z = globalPos.z();
  m_r = std::sqrt(m_x * m_x + m_y * m_y);
  m_varianceR = variance.x();
  m_varianceZ = variance.y();
}

template <typename external_spacepoint_t>
inline InternalSpacePoint<external_spacepoint_t>::InternalSpacePoint(
    const InternalSpacePoint<external_spacepoint_t>& sp)
    : m_sp(sp.sp()) {
  m_x = sp.m_x;
  m_y = sp.m_y;
  m_z = sp.m_z;
  m_r = sp.m_r;
  m_varianceR = sp.m_varianceR;
  m_varianceZ = sp.m_varianceZ;
}

}  // namespace Acts
