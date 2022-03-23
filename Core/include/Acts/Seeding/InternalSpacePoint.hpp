// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <array>
#include <cmath>
#include <limits>

namespace Acts {
template <typename SpacePoint>
class InternalSpacePoint {
  /////////////////////////////////////////////////////////////////////////////////
  // Public methods:
  /////////////////////////////////////////////////////////////////////////////////

 public:
  InternalSpacePoint() = delete;
  InternalSpacePoint(const SpacePoint& sp, const Acts::Vector3& globalPos,
                     const Acts::Vector2& offsetXY,
                     const Acts::Vector2& variance);

  InternalSpacePoint(const InternalSpacePoint<SpacePoint>& sp);
  ~InternalSpacePoint() = default;

  InternalSpacePoint<SpacePoint>& operator=(
      const InternalSpacePoint<SpacePoint>&) = delete;

  const float& x() const { return m_x; }
  const float& y() const { return m_y; }
  const float& z() const { return m_z; }
  const float& radius() const { return m_r; }
  float phi() const { return atan2f(m_y, m_x); }
  const float& varianceR() const { return m_varianceR; }
  const float& varianceZ() const { return m_varianceZ; }
  const float& deltaR() const { return m_deltaR; }
  const float& quality() const { return m_quality; }
  const float& cotTheta() const { return m_cotTheta; }
  void setDeltaR(float deltaR) { m_deltaR = deltaR; }
  void setCotTheta(float cotTheta) { m_cotTheta = cotTheta; }
  void setQuality(float quality) {
    if (quality >= m_quality) {
      m_quality = quality;
    }
  }
  const SpacePoint& sp() const { return m_sp; }

 protected:
  float m_x;          // x-coordinate in beam system coordinates
  float m_y;          // y-coordinate in beam system coordinates
  float m_z;          // z-coordinate in beam system coordinetes
  float m_r;          // radius       in beam system coordinates
  float m_varianceR;  //
  float m_varianceZ;  //
  float m_deltaR;     //
  float m_cotTheta = std::numeric_limits<
      double>::quiet_NaN();  // 1/tanTheta estimated from central+this space
                             // point. Its evaluation requires that the space
                             // point is a candidate for triplet search.
  float m_quality = -std::numeric_limits<
      double>::infinity();  // Quality score of the seed the space point is used
                            // for. Quality can be changed if the space point is
                            // used for a better quality seed.
  const SpacePoint& m_sp;   // external space point
};

/////////////////////////////////////////////////////////////////////////////////
// Inline methods
/////////////////////////////////////////////////////////////////////////////////

template <typename SpacePoint>
inline InternalSpacePoint<SpacePoint>::InternalSpacePoint(
    const SpacePoint& sp, const Acts::Vector3& globalPos,
    const Acts::Vector2& offsetXY, const Acts::Vector2& variance)
    : m_sp(sp) {
  m_x = globalPos.x() - offsetXY.x();
  m_y = globalPos.y() - offsetXY.y();
  m_z = globalPos.z();
  m_r = std::sqrt(m_x * m_x + m_y * m_y);
  m_varianceR = variance.x();
  m_varianceZ = variance.y();
}

/////////////////////////////////////////////////////////////////////////////////
// Copy constructor
/////////////////////////////////////////////////////////////////////////////////

template <typename SpacePoint>
inline InternalSpacePoint<SpacePoint>::InternalSpacePoint(
    const InternalSpacePoint<SpacePoint>& sp) = default;

}  // end of namespace Acts
