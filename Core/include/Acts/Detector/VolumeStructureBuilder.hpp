// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/interface/IExternalStructureBuilder.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <optional>

namespace Acts {
namespace Experimental {

/// This class provides the external detector volume structure, configured
/// either from:
///  - a volume extent
///  - from an array with volume bounds identification
///
class VolumeStructureBuilder : public IExternalStructureBuilder {
 public:
  /// Nexted configuration struct
  struct Config {
    /// A defined volume bound type
    VolumeBounds::BoundsType boundsType = VolumeBounds::BoundsType::eOther;
    /// The values (if already defined)
    std::vector<ActsScalar> boundValues = {};
    /// The optional extent to feed into the values
    std::optional<Extent> extent = std::nullopt;
  };

  /// Constructor
  ///
  /// @param cfg is the configuration struct
  /// @param logger logging instance for screen output
  VolumeStructureBuilder(const Config& cfg,
                         std::unique_ptr<const Logger> logger =
                             getDefaultLogger("VolumeStructureBuilder",
                                              Logging::INFO));

  /// The interface definition for internal structure creation
  ///
  /// @param gctx the geometry context at the creation of the internal structure
  ///
  /// @return a consistent set of detector volume internals
  ExternalStructure construct(const GeometryContext& gctx) const final;

 private:
  /// configuration object
  Config m_cfg;

  /// Private access method to the logger
  const Logger& logger() const { return *m_logger; }

  /// logging instance
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Experimental
}  // namespace Acts
