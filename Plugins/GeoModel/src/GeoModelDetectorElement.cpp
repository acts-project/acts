// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/GeoModel/GeoModelDetectorElement.hpp"

#include "Acts/Surfaces/Surface.hpp"

#include <utility>

#include <GeoModelKernel/GeoFullPhysVol.h>

using namespace Acts;

ActsPlugins::GeoModelDetectorElement::GeoModelDetectorElement(
    PVConstLink geoPhysVol, std::shared_ptr<Surface> surface,
    const Transform3& sfTransform, double thickness)
    : m_geoPhysVol(std::move(geoPhysVol)),
      m_surface(std::move(surface)),
      m_surfaceTransform(sfTransform),
      m_thickness(thickness) {
  if (m_surface) {
    attachSurface(std::move(m_surface));
  }
}

const Transform3& ActsPlugins::GeoModelDetectorElement::localToGlobalTransform(
    const GeometryContext& /*gctx*/) const {
  return m_surfaceTransform;
}

const Transform3& ActsPlugins::GeoModelDetectorElement::nominalTransform()
    const {
  return m_surfaceTransform;
}

const Surface& ActsPlugins::GeoModelDetectorElement::surface() const {
  return *m_surface;
}

Surface& ActsPlugins::GeoModelDetectorElement::surface() {
  return *m_surface;
}

double ActsPlugins::GeoModelDetectorElement::thickness() const {
  return m_thickness;
}

PVConstLink ActsPlugins::GeoModelDetectorElement::physicalVolume() const {
  return m_geoPhysVol;
}

const std::string& ActsPlugins::GeoModelDetectorElement::logVolName() const {
  return m_geoPhysVol->getLogVol()->getName();
}
