// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Volume.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"

#include <iostream>
#include <utility>

using namespace Acts::UnitLiterals;

Acts::Volume::Volume(const Transform3& transform,
                     std::shared_ptr<const VolumeBounds> volbounds)
    : GeometryObject(),
      m_transform(transform),
      m_itransform(m_transform.inverse()),
      m_center(m_transform.translation()),
      m_volumeBounds(std::move(volbounds)),
      m_orientedBoundingBox(m_volumeBounds->boundingBox(
          nullptr, {0.05_mm, 0.05_mm, 0.05_mm}, this)) {}

Acts::Volume::Volume(const Volume& vol, const Transform3& shift)
    : GeometryObject(),
      m_transform(shift * vol.m_transform),
      m_itransform(m_transform.inverse()),
      m_center(m_transform.translation()),
      m_volumeBounds(vol.m_volumeBounds),
      m_orientedBoundingBox(m_volumeBounds->boundingBox(
          nullptr, {0.05_mm, 0.05_mm, 0.05_mm}, this)) {}

Acts::Vector3 Acts::Volume::binningPosition(const GeometryContext& /*gctx*/,
                                            Acts::BinningValue bValue) const {
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

bool Acts::Volume::inside(const Acts::Vector3& gpos, double tol) const {
  Acts::Vector3 posInVolFrame((transform().inverse()) * gpos);
  return (volumeBounds()).inside(posInVolFrame, tol);
}

std::ostream& Acts::operator<<(std::ostream& sl, const Acts::Volume& vol) {
  sl << "Volume with " << vol.volumeBounds() << std::endl;
  return sl;
}

Acts::Volume::BoundingBox Acts::Volume::boundingBox(
    const Vector3& envelope) const {
  return m_volumeBounds->boundingBox(&m_transform, envelope, this);
}

const Acts::Volume::BoundingBox& Acts::Volume::orientedBoundingBox() const {
  return m_orientedBoundingBox;
}
