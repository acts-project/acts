// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
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
template <typename SpacePoint>
class InternalSpacePoint
{

  /////////////////////////////////////////////////////////////////////////////////
  // Public methods:
  /////////////////////////////////////////////////////////////////////////////////

public:
  InternalSpacePoint() = delete;
  InternalSpacePoint(const SpacePoint&     sp,
                     const Acts::Vector3D& globalPos,
                     const Acts::Vector2D& offsetXY,
                     const Acts::Vector2D& cov);

  InternalSpacePoint(const InternalSpacePoint<SpacePoint>& sp);
  ~InternalSpacePoint() = default;

  InternalSpacePoint<SpacePoint>&
  operator=(const InternalSpacePoint<SpacePoint>&);

  const float&
  x() const
  {
    return m_x;
  }
  const float&
  y() const
  {
    return m_y;
  }
  const float&
  z() const
  {
    return m_z;
  }
  const float&
  radius() const
  {
    return m_r;
  }
  float
  phi() const
  {
    return atan2f(m_y, m_x);
  }
  const float&
  covr() const
  {
    return m_covr;
  }
  const float&
  covz() const
  {
    return m_covz;
  }
  const SpacePoint&
  sp() const
  {
    return m_sp;
  }

protected:
  float             m_x;     // x-coordinate in beam system coordinates
  float             m_y;     // y-coordinate in beam system coordinates
  float             m_z;     // z-coordinate in beam system coordinetes
  float             m_r;     // radius       in beam system coordinates
  float             m_covr;  //
  float             m_covz;  //
  const SpacePoint& m_sp;    // external space point
};

/////////////////////////////////////////////////////////////////////////////////
// Inline methods
/////////////////////////////////////////////////////////////////////////////////

template <typename SpacePoint>
inline InternalSpacePoint<SpacePoint>::InternalSpacePoint(
    const SpacePoint&       sp,
    const Acts::Vector3D&   globalPos,
    const Acts::Vector2D&   offsetXY,
    const Acts::Vector2D&   cov): m_sp(sp)
{
  m_x    = globalPos.x() - offsetXY.x();
  m_y    = globalPos.y() - offsetXY.y();
  m_z    = globalPos.z();
  m_r    = std::sqrt(m_x * m_x + m_y * m_y);
  m_covr = cov.x();
  m_covz = cov.y();
}

/////////////////////////////////////////////////////////////////////////////////
// Copy constructor
/////////////////////////////////////////////////////////////////////////////////

template <typename SpacePoint>
inline InternalSpacePoint<SpacePoint>::InternalSpacePoint(
    const InternalSpacePoint<SpacePoint>& sp): m_sp(sp.sp())
{
  m_x    = sp.m_x;
  m_y    = sp.m_y;
  m_z    = sp.m_z;
  m_r    = sp.m_r;
  m_covr = sp.m_covr;
  m_covz = sp.m_covz;
}

}  // end of namespace Acts
