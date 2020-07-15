// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Volume.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Geometry/Volume.hpp"

#include <iostream>
#include <utility>

#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/Units.hpp"

using namespace Acts::UnitLiterals;

Acts::Volume::Volume()
    : GeometryObject(),
      m_transform(nullptr),
      m_center(s_origin),
      m_volumeBounds(nullptr),
      m_orientedBoundingBox(BoundingBox(this, {0, 0, 0}, {0, 0, 0})) {}

Acts::Volume::Volume(const std::shared_ptr<const Transform3D>& htrans,
                     std::shared_ptr<const VolumeBounds> volbounds)
    : GeometryObject(),
      m_transform(htrans),
      m_itransform(m_transform ? m_transform->inverse()
                               : Transform3D::Identity()),
      m_center(s_origin),
      m_volumeBounds(std::move(volbounds)),
      m_orientedBoundingBox(m_volumeBounds->boundingBox(
          nullptr, {0.05_mm, 0.05_mm, 0.05_mm}, this)) {
  if (htrans) {
    m_center = htrans->translation();
  }
}

Acts::Volume::Volume(const Volume& vol, const Transform3D* shift)
    : GeometryObject(),
      m_transform(vol.m_transform),
      m_itransform(m_transform ? m_transform->inverse()
                               : Transform3D::Identity()),
      m_center(s_origin),
      m_volumeBounds(vol.m_volumeBounds),
      m_orientedBoundingBox(m_volumeBounds->boundingBox(
          nullptr, {0.05_mm, 0.05_mm, 0.05_mm}, this)) {
  // apply the shift if it exists
  if (shift != nullptr) {
    m_transform = std::make_shared<const Transform3D>(transform() * (*shift));
    // reset inverse
    m_itransform = m_transform->inverse();
  }
  // now set the center
  m_center = transform().translation();
}

Acts::Volume::~Volume() = default;

const Acts::Vector3D Acts::Volume::binningPosition(
    const GeometryContext& /*gctx*/, Acts::BinningValue bValue) const {
  // for most of the binning types it is actually the center,
  // just for R-binning types the
  if (bValue == Acts::binR || bValue == Acts::binRPhi) {
    // the binning Position for R-type may have an offset
    return (center() + m_volumeBounds->binningOffset(bValue));
  }
  // return the center
  return center();
}

// assignment operator
Acts::Volume& Acts::Volume::operator=(const Acts::Volume& vol) {
  if (this != &vol) {
    m_transform = vol.m_transform;
    m_center = vol.m_center;
    m_volumeBounds = vol.m_volumeBounds;
  }
  return *this;
}

bool Acts::Volume::inside(const Acts::Vector3D& gpos, double tol) const {
  if (!m_transform) {
    return (volumeBounds()).inside(gpos, tol);
  }
  Acts::Vector3D posInVolFrame((transform().inverse()) * gpos);
  return (volumeBounds()).inside(posInVolFrame, tol);
}

std::ostream& Acts::operator<<(std::ostream& sl, const Acts::Volume& vol) {
  sl << "Volume with " << vol.volumeBounds() << std::endl;
  return sl;
}

Acts::Volume::BoundingBox Acts::Volume::boundingBox(
    const Vector3D& envelope) const {
  return m_volumeBounds->boundingBox(m_transform.get(), envelope, this);
}

const Acts::Volume::BoundingBox& Acts::Volume::orientedBoundingBox() const {
  return m_orientedBoundingBox;
}
