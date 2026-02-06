// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Volume.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"

#include <iostream>
#include <utility>

namespace Acts {

using namespace UnitLiterals;

Volume::Volume(const Transform3& transform,
               std::shared_ptr<VolumeBounds> volbounds)
    : GeometryObject(),
      m_transform(transform),
      m_itransform(m_transform.inverse()),
      m_center(m_transform.translation()),
      m_volumeBounds(std::move(volbounds)) {}

Volume::Volume(const Volume& vol, const Transform3& shift)
    : Volume(vol.shifted(shift)) {}

Volume Volume::shifted(const Transform3& shift) const {
  return Volume(shift * m_transform, m_volumeBounds);
}

Vector3 Volume::referencePosition(const GeometryContext& gctx,
                                  AxisDirection aDir) const {
  // for most of the binning types it is actually the center,
  // just for R-binning types the
  if (aDir == AxisDirection::AxisR || aDir == AxisDirection::AxisRPhi) {
    // the binning Position for R-type may have an offset
    return (center(gctx) + m_volumeBounds->referenceOffset(aDir));
  }
  // return the center
  return center(gctx);
}

// assignment operator
Volume& Volume::operator=(const Volume& vol) {
  if (this != &vol) {
    m_transform = vol.m_transform;
    m_center = vol.m_center;
    m_volumeBounds = vol.m_volumeBounds;
  }
  return *this;
}

bool Volume::inside(const GeometryContext& gctx, const Vector3& gpos,
                    double tol) const {
  Vector3 posInVolFrame = globalToLocalTransform(gctx) * gpos;
  return volumeBounds().inside(posInVolFrame, tol);
}
bool Volume::inside(const Vector3& gpos, double tol) const {
  ACTS_PUSH_IGNORE_DEPRECATED()
  return volumeBounds().inside(itransform() * gpos, tol);
  ACTS_POP_IGNORE_DEPRECATED()
}

std::ostream& operator<<(std::ostream& sl, const Volume& vol) {
  sl << "Volume with " << vol.volumeBounds() << std::endl;
  return sl;
}

Volume::BoundingBox Volume::boundingBox(const Vector3& envelope) const {
  return m_volumeBounds->boundingBox(&m_transform, envelope, this);
}

Volume::BoundingBox Volume::orientedBoundingBox() const {
  return m_volumeBounds->boundingBox(nullptr, {0.05_mm, 0.05_mm, 0.05_mm},
                                     this);
}

void Volume::assignVolumeBounds(std::shared_ptr<VolumeBounds> volbounds) {
  m_volumeBounds = std::move(volbounds);
}

void Volume::update(const GeometryContext& /*gctx*/,
                    std::shared_ptr<VolumeBounds> volbounds,
                    std::optional<Transform3> transform,
                    const Logger& /*logger*/) {
  if (volbounds) {
    m_volumeBounds = std::move(volbounds);
  }
  if (transform.has_value()) {
    setTransform(*transform);
  }
}

const Transform3& Volume::localToGlobalTransform(
    const GeometryContext& /*gctx*/) const {
  return m_transform;
}
const Transform3& Volume::globalToLocalTransform(
    const GeometryContext& /*gctx*/) const {
  return m_itransform;
}
const Transform3& Volume::transform() const {
  return m_transform;
}

const Transform3& Volume::itransform() const {
  return m_itransform;
}

const Vector3& Volume::center(const GeometryContext& /*gctx*/) const {
  return m_center;
}

const Vector3& Volume::center() const {
  return m_center;
}

const VolumeBounds& Volume::volumeBounds() const {
  return *m_volumeBounds;
}

VolumeBounds& Volume::volumeBounds() {
  return *m_volumeBounds;
}

std::shared_ptr<const VolumeBounds> Volume::volumeBoundsPtr() const {
  return m_volumeBounds;
}

std::shared_ptr<VolumeBounds> Volume::volumeBoundsPtr() {
  return m_volumeBounds;
}

void Volume::setTransform(const Transform3& transform) {
  m_transform = transform;
  m_itransform = m_transform.inverse();
  m_center = m_transform.translation();
}

bool Volume::operator==(const Volume& other) const {
  return (m_transform.matrix() == other.m_transform.matrix()) &&
         (*m_volumeBounds == *other.m_volumeBounds);
}

void Volume::visualize(IVisualization3D& helper, const GeometryContext& gctx,
                       const ViewConfig& viewConfig) const {
  auto bSurfaces =
      volumeBounds().orientedSurfaces(localToGlobalTransform(gctx));
  for (const auto& bs : bSurfaces) {
    bs.surface->visualize(helper, gctx, viewConfig);
  }
}

}  // namespace Acts
