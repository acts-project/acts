// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Traccc/TracccBackend.hpp"

#include <memory>

#include <traccc/bfield/magnetic_field.hpp>

namespace ActsExamples {

/// Magnetic field set up once for a traccc device backend
///
/// Holds an already constructed (covfie) host magnetic field and its copy on
/// the device. It is shared between all traccc device algorithms of a
/// sequence.
class TracccMagneticField {
 public:
  struct Config {
    /// The backend to set the field up for
    std::shared_ptr<const TracccBackend> backend;
    /// The host (covfie) magnetic field
    std::shared_ptr<const traccc::magnetic_field> field;
  };

  TracccMagneticField(const Config& cfg,
                      std::unique_ptr<const Acts::Logger> logger);

  const Config& config() const { return m_cfg; }

  const traccc::magnetic_field& hostField() const { return *m_cfg.field; }
  const traccc::magnetic_field& deviceField() const { return m_deviceField; }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  traccc::magnetic_field m_deviceField;
};

}  // namespace ActsExamples
