// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2019-07-17 Initial version
/// @author Moritz Kiehn <msmk@cern.ch>

#pragma once

#include <Acts/Utilities/Logger.hpp>
#include <memory>
#include <string>

#include "ACTFW/Framework/IService.hpp"

namespace FW {

/// A helper class for users to implement framework services.
///
/// This class provides default implementations for all interface methods and
/// and adds a default logger that can be used directly in subclasses.
/// Service implementations only need to implement the method that are
/// actually doing something.
class BareService : public IService {
 public:
  BareService(std::string name,
              Acts::Logging::Level level = Acts::Logging::INFO);

  /// The service name.
  std::string name() const final override;

  /// Default noop implementation for the start-of-run hook.
  void startRun() override;

  /// Default noop implementation for the per-event prepare hook.
  void prepare(AlgorithmContext& ctx) override;

 protected:
  const Acts::Logger& logger() const { return *m_logger; }

 private:
  std::string m_name;
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace FW
