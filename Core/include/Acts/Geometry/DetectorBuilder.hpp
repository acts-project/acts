// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ProtoDetector.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <memory>

namespace Acts {

namespace Experimental {
class Detector;

/// A Detector builder following a blueprint from a ProtoDetector
///
/// It takes some helper tools and a vector of surface objects,
/// together with a ProtoDetector description that is used to query a
/// KDTree for contained surfaces in structures defined by the proto
/// volume.
///
class DetectorBuilder {
 public:
  /// Nested Configuration for this TrackingGeometryBuilder
  struct Config {
    ProtoDetector protoDetector;
    /// Screen logging level for helper functions
    Logging::Level logLevel = Logging::INFO;
  };

  /// Constructor
  ///
  /// @param [in] cfg is the configuration struct for this builder
  /// @param [in] logger logging instance
  DetectorBuilder(const Config& cfg,
                  std::unique_ptr<const Logger> logger =
                      getDefaultLogger("DetectorBuilder", Logging::INFO));

  /// Destructor
  ~DetectorBuilder() = default;

  /// Detector building method
  ///
  /// @param gctx geometry context of that building call
  ///
  /// @return a shared pointer to the newly built detecots
  std::shared_ptr<Experimental::Detector> construct(
      const GeometryContext& gctx) const;

 private:
  /// Configuration member
  Config m_cfg;

  /// Private access method to the logger
  const Logger& logger() const { return *m_logger; }

  /// the logging instance
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Experimental
}  // namespace Acts
