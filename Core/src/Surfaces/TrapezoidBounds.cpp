// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"

#include <iomanip>
#include <iostream>

Acts::TrapezoidBounds::~TrapezoidBounds() = default;

Acts::SurfaceBounds::BoundsType Acts::TrapezoidBounds::type() const {
  return SurfaceBounds::eTrapezoid;
}

bool Acts::TrapezoidBounds::inside(const Acts::Vector2& lposition,
                                   const Acts::BoundaryCheck& bcheck) const {
  const double x = lposition[0];
  const double y = lposition[1];

  const double hlY = get(TrapezoidBounds::eHalfLengthY);
  const double hlXnY = get(TrapezoidBounds::eHalfLengthXnegY);
  const double hlXpY = get(TrapezoidBounds::eHalfLengthXposY);

  if (bcheck.type() == BoundaryCheck::Type::eAbsolute) {
    double tolX = bcheck.tolerance()[eBoundLoc0];
    double tolY = bcheck.tolerance()[eBoundLoc1];

    double distY = std::abs(y) - hlY;

    // check outside y range + tolerance
    if (distY > tolY) {
      return false;
    }

    // distance x to outer rectangle
    double distOuterX = std::abs(x) - std::max(hlXnY, hlXpY);

    // check outside x range + tolerance
    if (distOuterX > tolX) {
      return false;
    }

    // distance x to inner rectangle
    double distInnerX = std::abs(x) - std::min(hlXnY, hlXpY);

    // check inside x range, inside y range + tolerance (already checked above)
    if (distInnerX <= 0) {
      return true;
    }

    // check inside x range + tolerance, inside y range
    if ((distInnerX <= tolX) && (distY <= 0)) {
      return true;
    }
  }

  // at this stage, the point can only be in the triangles
  // run slow-ish polygon check
  std::array<Vector2, 4> v{
      Vector2{-hlXnY, -hlY}, {hlXnY, -hlY}, {hlXpY, hlY}, {-hlXpY, hlY}};
  return bcheck.isInside(lposition, v);
}

std::vector<Acts::Vector2> Acts::TrapezoidBounds::vertices(
    unsigned int /*lseg*/) const {
  double minhx = get(TrapezoidBounds::eHalfLengthXnegY);
  double maxhx = get(TrapezoidBounds::eHalfLengthXposY);
  double hy = get(TrapezoidBounds::eHalfLengthY);
  return {{-minhx, -hy}, {minhx, -hy}, {maxhx, hy}, {-maxhx, hy}};
}

const Acts::RectangleBounds& Acts::TrapezoidBounds::boundingBox() const {
  return m_boundingBox;
}

std::ostream& Acts::TrapezoidBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::TrapezoidBounds:  (halfXnegY, halfXposY, halfY) = "
     << "(" << get(eHalfLengthXnegY) << ", " << get(eHalfLengthXposY) << ", "
     << get(eHalfLengthY) << ")";
  sl << std::setprecision(-1);
  return sl;
}
