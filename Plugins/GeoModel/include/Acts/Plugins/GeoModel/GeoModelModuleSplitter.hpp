// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>
#include <string>
#include <vector>

namespace Acts {

/// This class allows to split modules with annulus shape along into several
/// parts with the help of a list of radii. This is used for building
/// an ITk geometry.
class GeoModelModuleSplitter {
 public:
  using DetectorElementPtr = std::shared_ptr<GeoModelDetectorElement>;

  /// Construct the module splitter
  /// @param splitPatterns A map of named split patterns
  /// @param tolerance The tolerance applied when matching the split
  /// patterns to the modules
  /// @param level The level of the internal logger
  GeoModelModuleSplitter(
      const std::map<std::string, std::vector<double>> &splitPatterns,
      double tolerance = 1.e-3,
      Acts::Logging::Level level = Acts::Logging::INFO)
      : m_splitPatterns(splitPatterns),
        m_tolerance(tolerance),
        m_logger(getDefaultLogger("GeoModelModuleSplitter", level)) {}

  /// Perform the split for a detector element
  /// @param detElement The detector element to split
  /// @param gctx The geometry context
  std::optional<std::vector<DetectorElementPtr>> split(
      DetectorElementPtr detElement, const GeometryContext &gctx) const;

 private:
  std::map<std::string, std::vector<double>> m_splitPatterns;
  double m_tolerance;
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace Acts
