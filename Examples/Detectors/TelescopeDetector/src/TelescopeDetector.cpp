// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "ActsExamples/TelescopeDetector/BuildTelescopeDetector.hpp"

#include <stdexcept>

namespace ActsExamples {

TelescopeDetector::TelescopeDetector(const Config& cfg)
    : Detector(Acts::getDefaultLogger("TelescopeDetector", cfg.logLevel)),
      m_cfg(cfg) {
  if (m_cfg.surfaceType > 1) {
    throw std::invalid_argument(
        "The surface type could either be 0 for plane surface or 1 for disc "
        "surface.");
  }
  if (m_cfg.rotDirection > 2) {
    throw std::invalid_argument("The axis value could only be 0, 1, or 2.");
  }
  // Check if the bounds values are valid
  if (m_cfg.surfaceType == 1 && m_cfg.bounds[0] >= m_cfg.bounds[1]) {
    throw std::invalid_argument(
        "The minR should be smaller than the maxR for disc surface bounds.");
  }

  if (m_cfg.positions.empty()) {
    throw std::invalid_argument("At least one surface position is required.");
  }

  if (m_cfg.positions.size() != m_cfg.stereos.size()) {
    throw std::invalid_argument(
        "The number of provided positions must match the number of "
        "provided stereo angles.");
  }

  m_nominalGeometryContext =
      Acts::GeometryContext::dangerouslyDefaultConstruct();

  if (m_cfg.gen3) {
    m_trackingGeometry = buildTelescopeDetectorGen3(
        m_nominalGeometryContext, m_detectorStore, m_cfg.positions,
        m_cfg.stereos, m_cfg.offsets, m_cfg.bounds, m_cfg.thickness,
        static_cast<TelescopeSurfaceType>(m_cfg.surfaceType),
        static_cast<Acts::AxisDirection>(m_cfg.rotDirection), m_cfg.envelope_x,
        m_cfg.envelope_y, m_cfg.envelope_z, logger());
  } else {
    m_trackingGeometry = buildTelescopeDetector(
        m_nominalGeometryContext, m_detectorStore, m_cfg.positions,
        m_cfg.stereos, m_cfg.offsets, m_cfg.bounds, m_cfg.thickness,
        static_cast<TelescopeSurfaceType>(m_cfg.surfaceType),
        static_cast<Acts::AxisDirection>(m_cfg.rotDirection));
  }
}

TelescopeDetector::TelescopeDetector(const Config& cfg, NoBuildTag /*unused*/)
    : Detector(Acts::getDefaultLogger("TelescopeDetector", cfg.logLevel)),
      m_cfg(cfg) {
  if (m_cfg.surfaceType > 1) {
    throw std::invalid_argument(
        "The surface type could either be 0 for plane surface or 1 for disc "
        "surface.");
  }
  if (m_cfg.rotDirection > 2) {
    throw std::invalid_argument("The axis value could only be 0, 1, or 2.");
  }
  // Check if the bounds values are valid
  if (m_cfg.surfaceType == 1 && m_cfg.bounds[0] >= m_cfg.bounds[1]) {
    throw std::invalid_argument(
        "The minR should be smaller than the maxR for disc surface bounds.");
  }

  if (m_cfg.positions.empty()) {
    throw std::invalid_argument("At least one surface position is required.");
  }

  if (m_cfg.positions.size() != m_cfg.stereos.size()) {
    throw std::invalid_argument(
        "The number of provided positions must match the number of "
        "provided stereo angles.");
  }
}

std::shared_ptr<Acts::PlaneSurface> TelescopeDetector::getReferenceSurface(
    const Acts::Vector3& position, double halfX, double halfY) const {
  using namespace Acts::UnitLiterals;

  Acts::RotationMatrix3 rotation = Acts::RotationMatrix3::Identity();
  if (m_cfg.rotDirection == 0) {  // 0 == Acts::AxisDirection::AxisX
    Acts::AngleAxis3 rot{90._degree, Acts::Vector3::UnitY()};
    rotation = rot.toRotationMatrix();
  } else if (m_cfg.rotDirection == 1) {  // 1 == Acts::AxisDirection::AxisY
    Acts::AngleAxis3 rot{90._degree, Acts::Vector3::UnitX()};
    rotation = rot.toRotationMatrix();
  }

  Acts::Transform3 transform = Acts::Transform3::Identity();
  transform.linear() = rotation;
  transform.translation() = position;

  auto bounds = std::make_shared<Acts::RectangleBounds>(halfX, halfY);
  return Acts::Surface::makeShared<Acts::PlaneSurface>(transform, bounds);
}

}  // namespace ActsExamples
