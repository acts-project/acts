// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/interface/IDetectorComponentBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>
#include <vector>

namespace Acts {
namespace Experimental {

class IRootVolumeFinderBuilder;
class IGeometryIdGenerator;

/// @brief A dedicated container builder for cylindrical detector containers
///
/// It relies on the detailed implementation of the CylindricalDetectorHelper
/// and allows for DetectorVolume attachment in z/r/phi, such as wrapping
/// of bevelled cylinder objects in z/r
///
/// @note the builder expects a fully consistent set of sub volume builders
/// that will be executed in a chain
///
/// @note allowed BinningValue(s) for the cylindrical container builder are
/// {binZ}, {binR}, {binPhi}, {binZ, binR}, whereas the last option indicates
/// a wrapping setup.
class CylindricalContainerBuilder : public IDetectorComponentBuilder {
 public:
  /// Nested configuration object
  struct Config {
    /// The configured volume builders
    std::vector<std::shared_ptr<const IDetectorComponentBuilder>> builders = {};
    /// Binning prescription of attachment
    std::vector<BinningValue> binning = {};
    /// The root volume finder
    std::shared_ptr<IRootVolumeFinderBuilder> rootVolumeFinderBuilder = nullptr;
    /// The geometry id generator
    std::shared_ptr<const IGeometryIdGenerator> geoIdGenerator = nullptr;
    /// An eventual reverse geometry id generation
    bool geoIdReverseGen = false;
    /// Auxiliary information, mainly for screen output
    std::string auxiliary = "";
  };

  /// Constructor with configuration arguments
  ///
  /// @param cfg is the configuration struct
  /// @param logger logging instance for screen output
  CylindricalContainerBuilder(
      const Config& cfg,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("CylindricalContainerBuilder", Logging::INFO));

  /// The final implementation of the cylindrical container builder
  ///
  /// @param gctx The geometry context for this call
  ///
  /// @return an outgoing detector component
  DetectorComponent construct(const GeometryContext& gctx) const final;

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
