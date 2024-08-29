// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"

#include "Acts/Surfaces/Surface.hpp"

#include <utility>

Acts::GeoModelDetectorElement::GeoModelDetectorElement(
    PVConstLink geoPhysVol, std::shared_ptr<Surface> surface,
    const Transform3& sfTransform, ActsScalar thickness)
    : m_geoPhysVol(std::move(geoPhysVol)),
      m_surface(std::move(surface)),
      m_surfaceTransform(sfTransform),
      m_thickness(thickness) {}

const Acts::Transform3& Acts::GeoModelDetectorElement::transform(
    const GeometryContext& /*gctx*/) const {
  return m_surfaceTransform;
}

const Acts::Surface& Acts::GeoModelDetectorElement::surface() const {
  return *m_surface;
}

Acts::Surface& Acts::GeoModelDetectorElement::surface() {
  return *m_surface;
}

Acts::ActsScalar Acts::GeoModelDetectorElement::thickness() const {
  return m_thickness;
}

PVConstLink Acts::GeoModelDetectorElement::physicalVolume() const {
  return m_geoPhysVol;
}
