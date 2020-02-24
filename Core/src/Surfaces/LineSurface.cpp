// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/LineSurface.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

#include "Acts/Utilities/ThrowAssert.hpp"

Acts::LineSurface::LineSurface(std::shared_ptr<const Transform3D> htrans,
                               double radius, double halez)
    : GeometryObject(),
      Surface(std::move(htrans)),
      m_bounds(std::make_shared<const LineBounds>(radius, halez)) {}

Acts::LineSurface::LineSurface(std::shared_ptr<const Transform3D> htrans,
                               std::shared_ptr<const LineBounds> lbounds)
    : GeometryObject(),
      Surface(std::move(htrans)),
      m_bounds(std::move(lbounds)) {}

Acts::LineSurface::LineSurface(const std::shared_ptr<const LineBounds>& lbounds,
                               const DetectorElementBase& detelement)
    : GeometryObject(), Surface(detelement), m_bounds(lbounds) {
  throw_assert(lbounds, "LineBounds must not be nullptr");
}

Acts::LineSurface::LineSurface(const LineSurface& other)
    : GeometryObject(), Surface(other), m_bounds(other.m_bounds) {}

Acts::LineSurface::LineSurface(const GeometryContext& gctx,
                               const LineSurface& other,
                               const Transform3D& transf)
    : GeometryObject(),
      Surface(gctx, other, transf),
      m_bounds(other.m_bounds) {}

Acts::LineSurface& Acts::LineSurface::operator=(const LineSurface& other) {
  if (this != &other) {
    Surface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}
