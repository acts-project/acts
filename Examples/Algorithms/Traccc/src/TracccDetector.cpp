// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/TracccDetector.hpp"

#include <stdexcept>

namespace ActsExamples {

TracccDetector::TracccDetector(const Config& cfg,
                               std::unique_ptr<const Acts::Logger> _logger)
    : m_cfg(cfg), m_logger(std::move(_logger)) {
  if (m_logger == nullptr) {
    throw std::invalid_argument("TracccDetector: logger is not set");
  }
  if (m_cfg.backend == nullptr) {
    throw std::invalid_argument("TracccDetector: backend is not set");
  }
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("TracccDetector: trackingGeometry is not set");
  }
  if (m_cfg.detrayDetector == nullptr) {
    throw std::invalid_argument("TracccDetector: detrayDetector is not set");
  }
  if (m_cfg.surfaceIds == nullptr) {
    throw std::invalid_argument("TracccDetector: surfaceIds is not set");
  }

  const TracccBackend& backend = *m_cfg.backend;

  // Copy the detector to the device
  m_deviceDetector = traccc::buffer_from_host_detector(
      *m_cfg.detrayDetector, backend.memoryResource().main, backend.copy());
  backend.synchronize();
  ACTS_INFO("Copied detray detector with " << m_cfg.surfaceIds->size()
                                           << " sensitive surfaces to the "
                                           << backend.name() << " device");
}

}  // namespace ActsExamples
