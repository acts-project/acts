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
Volume::Volume(const Transform3& transform,
               std::shared_ptr<VolumeBounds> volbounds) noexcept
    : GeometryObject(),
      m_transform{std::make_unique<Transform3>(transform)},
      m_itransform{std::make_unique<Transform3>(transform.inverse())},
      m_center{transform.translation()},
      m_volumeBounds(std::move(volbounds)) {}

Volume Volume::shifted(const GeometryContext& gctx,
                       const Transform3& shift) const {
  return Volume(shift * localToGlobalTransform(gctx), m_volumeBounds);
}
Volume::Volume(VolumePlacementBase& positioner,
               std::shared_ptr<VolumeBounds> volbounds) noexcept
    : GeometryObject{},
      m_volumeBounds{std::move(volbounds)},
      m_placement{&positioner} {}

Volume::Volume(const Volume& vol, const Transform3& shift)
    : Volume{shift * (vol.m_transform ? (*vol.m_transform)
                                      : Transform3::Identity()),
             vol.m_volumeBounds} {}

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

bool Volume::isAlignable() const {
  return m_placement != nullptr;
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
  return m_volumeBounds->boundingBox(m_transform.get(), envelope, this);
}

Volume::BoundingBox Volume::orientedBoundingBox() const {
  using namespace UnitLiterals;
  return m_volumeBounds->boundingBox(nullptr, {0.05_mm, 0.05_mm, 0.05_mm},
                                     this);
}

void Volume::assignVolumeBounds(std::shared_ptr<VolumeBounds> volbounds) {
  assert(volbounds != nullptr);
  // If the volume is instantiated with a placement, the bounds can be updated
  // as long as the portals have not been made. Or the bounds are equivalent
  // with the current bounds
  if (isAlignable() && volumePlacement()->nPortalPlacements() > 0ul &&
      (*m_volumeBounds) != (*volbounds)) {
    throw std::runtime_error(
        "assignVolumeBounds() - Bounds cannot be overwritten if the associated "
        "VolumePlacement has instantiated boundary surfaces");
  }

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
    const GeometryContext& gctx) const {
  if (isAlignable()) {
    return volumePlacement()->localToGlobalTransform(gctx);
  }
  assert(m_transform != nullptr);
  return (*m_transform);
}
const Transform3& Volume::globalToLocalTransform(
    const GeometryContext& gctx) const {
  if (isAlignable()) {
    return volumePlacement()->globalToLocalTransform(gctx);
  }
  assert(m_itransform != nullptr);
  return (*m_itransform);
}
const Transform3& Volume::transform() const {
  assert(m_transform != nullptr);
  return (*m_transform);
}

const Transform3& Volume::itransform() const {
  assert(m_itransform != nullptr);
  return (*m_itransform);
}

Vector3 Volume::center(const GeometryContext& gctx) const {
  return localToGlobalTransform(gctx).translation();
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

VolumePlacementBase* Volume::volumePlacement() {
  return m_placement;
}
const VolumePlacementBase* Volume::volumePlacement() const {
  return m_placement;
}

void Volume::setTransform(const Transform3& transform) {
  if (isAlignable()) {
    throw std::runtime_error(
        "setTransform() - Transforms of externally aligned volumes cannot "
        "be overwritten");
  }
  m_transform = std::make_unique<Transform3>(transform);
  m_itransform = std::make_unique<Transform3>(transform.inverse());
  m_center = transform.translation();
}

bool Volume::operator==(const Volume& other) const {
  return ((m_transform != nullptr && other.m_transform != nullptr &&
           m_transform->matrix() == other.m_transform->matrix()) ||
          (volumePlacement() == other.volumePlacement() && isAlignable())) &&
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
