// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cmath>
#include <functional>
#include <optional>

namespace Acts {

template <typename SpacePoint>
class InternalSpacePoint {
  /////////////////////////////////////////////////////////////////////////////////
  // Public methods:
  /////////////////////////////////////////////////////////////////////////////////

 public:
  InternalSpacePoint() = delete;
  InternalSpacePoint(std::size_t index, const SpacePoint& sp,
                     const Acts::Vector3& globalPos,
                     const Acts::Vector2& offsetXY,
                     const Acts::Vector2& variance,
                     std::optional<float> globalTime);

  InternalSpacePoint(const InternalSpacePoint<SpacePoint>& sp);
  ~InternalSpacePoint() = default;

  InternalSpacePoint<SpacePoint>& operator=(
      const InternalSpacePoint<SpacePoint>&) = delete;

  std::size_t index() const { return m_index; }
  float x() const { return m_x; }
  float y() const { return m_y; }
  float z() const { return m_z; }
  std::optional<float> t() const { return m_t; }
  float radius() const { return m_r; }
  float phi() const { return m_phi; }
  float varianceR() const { return m_varianceR; }
  float varianceZ() const { return m_varianceZ; }
  const SpacePoint& sp() const { return m_sp; }

 protected:
  std::size_t m_index;
  float m_x;                 // x-coordinate in beam system coordinates
  float m_y;                 // y-coordinate in beam system coordinates
  float m_z;                 // z-coordinate in beam system coordinetes
  float m_r;                 // radius       in beam system coordinates
  float m_phi;               //
  float m_varianceR;         //
  float m_varianceZ;         //
  std::optional<float> m_t;  // time
  std::reference_wrapper<const SpacePoint> m_sp;  // external space point
};

/////////////////////////////////////////////////////////////////////////////////
// Inline methods
/////////////////////////////////////////////////////////////////////////////////

template <typename SpacePoint>
inline InternalSpacePoint<SpacePoint>::InternalSpacePoint(
    std::size_t index, const SpacePoint& sp, const Acts::Vector3& globalPos,
    const Acts::Vector2& offsetXY, const Acts::Vector2& variance,
    std::optional<float> globalTime)
    : m_index(index),
      m_x(globalPos.x() - offsetXY.x()),
      m_y(globalPos.y() - offsetXY.y()),
      m_z(globalPos.z()),
      m_r(std::hypot(m_x, m_y)),
      m_phi(std::atan2(m_y, m_x)),
      m_varianceR(variance.x()),
      m_varianceZ(variance.y()),
      m_t(globalTime),
      m_sp(sp) {}

template <typename SpacePoint>
inline bool operator<(const InternalSpacePoint<SpacePoint>& lhs,
                      const InternalSpacePoint<SpacePoint>& rhs) {
  // TODO would it be sufficient to check just the index under the assumption
  //   that the same measurement index always produces the same space point?
  // no need to check r since it is fully defined by x/y
  // Can be used for sorting the space points
  if ((lhs.x() == rhs.x()) and (lhs.y() == rhs.y()) and (lhs.z() == rhs.z()) and
      (lhs.varianceR() == rhs.varianceR()) and
      (lhs.varianceZ() == rhs.varianceZ())) {
    return false;
  }
  if (lhs.x() != rhs.x()) {
    return (lhs.x() < rhs.x());
  }
  if (lhs.y() != rhs.y()) {
    return (lhs.y() < rhs.y());
  }
  if (lhs.z() != rhs.z()) {
    return (lhs.z() < rhs.z());
  }
  return false;
}

/////////////////////////////////////////////////////////////////////////////////
// Copy constructor
/////////////////////////////////////////////////////////////////////////////////

template <typename SpacePoint>
inline InternalSpacePoint<SpacePoint>::InternalSpacePoint(
    const InternalSpacePoint<SpacePoint>& sp) = default;

}  // namespace Acts
