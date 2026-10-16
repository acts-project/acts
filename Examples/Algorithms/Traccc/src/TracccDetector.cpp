// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/TracccDetector.hpp"

#include <stdexcept>

#include <traccc/bfield/construct_const_bfield.hpp>
#include <traccc/io/read_magnetic_field.hpp>

namespace ActsExamples {

TracccDetector::TracccDetector(const Config& cfg,
                               std::unique_ptr<const Acts::Logger> _logger)
    : m_cfg(cfg),
      m_logger(_logger != nullptr
                   ? std::move(_logger)
                   : Acts::getDefaultLogger("TracccDetector",
                                            Acts::Logging::INFO)) {
  if (m_cfg.backend == nullptr) {
    throw std::invalid_argument("TracccDetector: backend is not set");
  }
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("TracccDetector: trackingGeometry is not set");
  }
  if (m_cfg.detrayDetector == nullptr) {
    throw std::invalid_argument("TracccDetector: detrayDetector is not set");
  }

  const TracccBackend& backend = *m_cfg.backend;
  const DetectorType& detector = *m_cfg.detrayDetector;

  // detray <-> ACTS identifier maps: the payload converter stores the ACTS
  // geometry identifier of every surface as its detray source link
  for (const auto& surface : detector.surfaces()) {
    if (!surface.is_sensitive()) {
      continue;
    }
    const Acts::GeometryIdentifier actsId(surface.source);
    const std::uint64_t detrayId = surface.identifier().value();
    m_detrayToActs.emplace(detrayId, actsId);
    m_actsToDetray.emplace(actsId, detrayId);
  }
  ACTS_INFO("Mapped " << m_detrayToActs.size()
                      << " sensitive surfaces between detray and ACTS");

  // Copy the detector to the device
  m_deviceDetector.set<DetectorType>(detray::get_buffer(
      detector, backend.memoryResource().main, backend.copy()));
  backend.synchronize();
  ACTS_INFO("Copied detray detector with "
            << detector.volumes().size() << " volumes and "
            << detector.surfaces().size() << " surfaces to the "
            << backend.name() << " device");

  // Magnetic field
  if (!m_cfg.bFieldFile.empty()) {
    ACTS_INFO("Reading magnetic field from " << m_cfg.bFieldFile);
    traccc::io::read_magnetic_field(m_hostField, m_cfg.bFieldFile);
  } else {
    ACTS_INFO("Using constant magnetic field "
              << m_cfg.constantBField.transpose());
    m_hostField = traccc::construct_const_bfield(
        static_cast<traccc::scalar>(m_cfg.constantBField.x()),
        static_cast<traccc::scalar>(m_cfg.constantBField.y()),
        static_cast<traccc::scalar>(m_cfg.constantBField.z()));
  }
  m_deviceField = backend.makeMagneticField(m_hostField);
}

std::optional<Acts::GeometryIdentifier> TracccDetector::actsIdentifier(
    std::uint64_t detrayId) const {
  auto it = m_detrayToActs.find(detrayId);
  if (it == m_detrayToActs.end()) {
    return std::nullopt;
  }
  return it->second;
}

std::optional<std::uint64_t> TracccDetector::detrayIdentifier(
    Acts::GeometryIdentifier actsId) const {
  auto it = m_actsToDetray.find(actsId);
  if (it == m_actsToDetray.end()) {
    return std::nullopt;
  }
  return it->second;
}

}  // namespace ActsExamples
