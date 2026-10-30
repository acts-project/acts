// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/TracccMagneticField.hpp"

#include <stdexcept>

namespace ActsExamples {

TracccMagneticField::TracccMagneticField(
    const Config& cfg, std::unique_ptr<const Acts::Logger> _logger)
    : m_cfg(cfg), m_logger(std::move(_logger)) {
  if (m_logger == nullptr) {
    throw std::invalid_argument("TracccMagneticField: logger is not set");
  }
  if (m_cfg.backend == nullptr) {
    throw std::invalid_argument("TracccMagneticField: backend is not set");
  }
  if (m_cfg.field == nullptr) {
    throw std::invalid_argument("TracccMagneticField: field is not set");
  }

  m_deviceField = m_cfg.backend->makeMagneticField(*m_cfg.field);
  ACTS_INFO("Copied magnetic field to the " << m_cfg.backend->name()
                                            << " device");
}

}  // namespace ActsExamples
