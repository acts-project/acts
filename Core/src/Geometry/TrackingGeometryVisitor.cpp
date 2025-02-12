// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometryVisitor.hpp"

namespace Acts {

TrackingGeometryVisitor::~TrackingGeometryVisitor() = default;

void TrackingGeometryVisitor::visitVolume(const TrackingVolume& /*volume*/) {
  // Default implementation is a no-op
}

void TrackingGeometryVisitor::visitPortal(const Portal& /*portal*/) {
  // Default implementation is a no-op
}

void TrackingGeometryVisitor::visitSurface(const Surface& /*surface*/) {
  // Default implementation is a no-op
}

void TrackingGeometryVisitor::visitLayer(const Layer& /*layer*/) {
  // Default implementation is a no-op
}

void TrackingGeometryVisitor::visitBoundarySurface(
    const BoundarySurfaceT<TrackingVolume>& /*boundary*/) {
  // Default implementation is a no-op
}

void TrackingGeometryMutableVisitor::visitVolume(TrackingVolume& /*volume*/) {
  // Default implementation is a no-op
}

void TrackingGeometryMutableVisitor::visitPortal(Portal& /*portal*/) {
  // Default implementation is a no-op
}

void TrackingGeometryMutableVisitor::visitSurface(Surface& /*surface*/) {
  // Default implementation is a no-op
}

void TrackingGeometryMutableVisitor::visitBoundarySurface(
    BoundarySurfaceT<TrackingVolume>& /*boundary*/) {
  // Default implementation is a no-op
}

void TrackingGeometryMutableVisitor::visitLayer(Layer& /*layer*/) {
  // Default implementation is a no-op
}

TrackingGeometryLambdaVisitor::TrackingGeometryLambdaVisitor(Config&& config)
    : m_config(std::move(config)) {}

void TrackingGeometryLambdaVisitor::visitVolume(const TrackingVolume& volume) {
  if (m_config.volume) {
    m_config.volume(volume);
  }
}

void TrackingGeometryLambdaVisitor::visitPortal(const Portal& portal) {
  if (m_config.portal) {
    m_config.portal(portal);
  }
}

void TrackingGeometryLambdaVisitor::visitSurface(const Surface& surface) {
  if (m_config.surface) {
    m_config.surface(surface);
  }
}

void TrackingGeometryLambdaVisitor::visitLayer(const Layer& layer) {
  if (m_config.layer) {
    m_config.layer(layer);
  }
}

void TrackingGeometryLambdaVisitor::visitBoundarySurface(
    const BoundarySurfaceT<TrackingVolume>& boundary) {
  if (m_config.boundary) {
    m_config.boundary(boundary);
  }
}

TrackingGeometryMutableVisitor::~TrackingGeometryMutableVisitor() = default;

TrackingGeometryLambdaMutableVisitor::TrackingGeometryLambdaMutableVisitor(
    Config&& config)
    : m_config(std::move(config)) {}

void TrackingGeometryLambdaMutableVisitor::visitVolume(TrackingVolume& volume) {
  if (m_config.volume) {
    m_config.volume(volume);
  }
}

void TrackingGeometryLambdaMutableVisitor::visitPortal(Portal& portal) {
  if (m_config.portal) {
    m_config.portal(portal);
  }
}

void TrackingGeometryLambdaMutableVisitor::visitSurface(Surface& surface) {
  if (m_config.surface) {
    m_config.surface(surface);
  }
}

void TrackingGeometryLambdaMutableVisitor::visitLayer(Layer& layer) {
  if (m_config.layer) {
    m_config.layer(layer);
  }
}

void TrackingGeometryLambdaMutableVisitor::visitBoundarySurface(
    BoundarySurfaceT<TrackingVolume>& boundary) {
  if (m_config.boundary) {
    m_config.boundary(boundary);
  }
}

}  // namespace Acts
