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

class ModuleByRadiusSplitter {
 public:
  /// Construct with a map { design-name: list of radii }
  ModuleByRadiusSplitter(
      const std::map<std::string, std::vector<double>> &designs,
      double tolerance = 1.e-3,
      Acts::Logging::Level level = Acts::Logging::INFO)
      : m_designs(designs),
        m_tolerance(tolerance),
        m_logger(getDefaultLogger("ModuleByRadiusSplitter", level)) {}

  std::optional<std::vector<std::pair<
      std::shared_ptr<Surface>, std::shared_ptr<GeoModelDetectorElement>>>>
  split(std::shared_ptr<Surface> surface,
        std::shared_ptr<GeoModelDetectorElement> detElement,
        const GeometryContext &gctx) const;

 private:
  std::map<std::string, std::vector<double>> m_designs;
  double m_tolerance;
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace Acts
