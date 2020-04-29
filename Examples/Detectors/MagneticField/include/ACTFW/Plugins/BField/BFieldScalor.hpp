// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

#include "ACTFW/Framework/AlgorithmContext.hpp"
#include "ACTFW/Framework/IContextDecorator.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace FW {

namespace BField {

/// A mockup service that rotates a
/// cylindrical geometry
class BFieldScalor : public IContextDecorator {
 public:
  /// @brief nested configuration struct
  struct Config {
    /// Incremental scaling
    double scalor = 1.2;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param logger The logging framework
  BFieldScalor(const Config& cfg,
               std::unique_ptr<const Acts::Logger> logger =
                   Acts::getDefaultLogger("BFieldScalor", Acts::Logging::INFO));

  /// Virtual destructor
  virtual ~BFieldScalor() = default;

  /// @brief decorates (adds, modifies) the AlgorithmContext
  /// with a geometric rotation per event
  ///
  /// @note If decorators depend on each other, they have to be
  /// added in order.
  ///
  /// @param context the bare (or at least non-const) Event context
  ProcessCode decorate(AlgorithmContext& context) final override;

  /// @brief decorator name() for screen output
  const std::string& name() const final override { return m_name; }

 private:
  Config m_cfg;                                  ///< the configuration class
  std::unique_ptr<const Acts::Logger> m_logger;  ///!< the logging instance
  std::string m_name = "BFieldScalor";

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }
};
}  // namespace BField

}  // namespace FW
