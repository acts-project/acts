// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/Blueprint.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>
#include <tuple>

namespace dd4hep {
class DetElement;
}

namespace Acts {

namespace Experimental {

/// @brief This class allows to generate layer structure builders for dd4hep sub detectors
/// It performs an intermediate step by taking dd4hep::DetElemnent objects that
/// describe a detector structure and prepares the translation of the detector
/// element and eventual passive surfaces.
///
/// It would also build passive support structures if configured to do so.
///
class DD4hepDetectorStructure {
 public:
  /// @brief  Nested options struct for the building
  struct Options {
    /// The log level of the tools
    Logging::Level logLevel = Logging::INFO;
    /// Boolean to flag if blueprint should be written out
    bool blueprintDot = false;
  };

  /// Constructor with from file name
  ///
  /// @param logger is the screen output logger
  ///
  /// @note this needs to be provided
  DD4hepDetectorStructure(std::unique_ptr<const Logger> logger =
                              getDefaultLogger("DD4hepLayerStructure",
                                               Acts::Logging::INFO));

  DD4hepDetectorStructure() = delete;

  /// @brief This method builds a DD4hep detector
  ///
  /// @param gctx the geometry context
  /// @param dd4hepElement is the dd4hep::DetElement of the world volume
  /// @param options is the options struct
  ///
  /// @return a detector, a blue node and a store
  std::tuple<std::shared_ptr<const Detector>, DD4hepDetectorElement::Store>
  construct(const GeometryContext& gctx,
            const dd4hep::DetElement& dd4hepElement,
            const Options& options) const;

 private:
  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to the logger
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Experimental
}  // namespace Acts
