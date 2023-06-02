// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/interface/IDetectorBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <vector>

namespace Acts {

class Surface;

namespace Experimental {

class Detector;

struct NoSurfaces {
  std::vector<std::shared_ptr<Surface>> operator()() const { return {}; }
};

/// @brief A dedicated container builder for detectors based on a grid structure
/// It supports the constructruction of:
///
/// - a cylindrical grid in z/r/phi (and subsets)
/// - a cartesian grid in x/y/z (and subsets)
///
/// @todo change to a component builder (after moving root volumes into the DetectorComponent)
///
class GridDetectorBuilder : public IDetectorBuilder {
 public:
  /// @brief The surface binning definition
  struct Binning {
    /// Define the binning of the surfaces
    BinningData data;
    /// An expansion for the filling (in bins)
    size_t expansion = 0u;
  };

  /// @brief Nested configuration struct
  struct Config {
    // Name of the detector builder
    std::string name = "";
    /// Connection point for a function to provide surfaces
    std::function<std::vector<std::shared_ptr<Surface>>()> surfaces =
        NoSurfaces{};
    /// Binning description
    std::vector<Binning> binning = {};
    /// The segments for the polyhedron reference generator
    unsigned int polyhedronSegements = 4;
    /// A gobal transform
    Transform3 transform = Transform3::Identity();
    /// Auxiliary information
    std::string auxillary = "";
  };

  /// Constructor
  ///
  /// @param cfg is the configuration struct
  /// @param logger logging instance for screen output
  GridDetectorBuilder(const Config& cfg,
                      std::unique_ptr<const Logger> logger = getDefaultLogger(
                          "GridDetectorBuilder", Logging::INFO));

  /// The interface definition for internal structure creation
  ///
  /// @param gctx the geometry context at the creation of the internal structure
  ///
  /// @return a consistent set of detector volume internals
  std::shared_ptr<Detector> construct(const GeometryContext& gctx) const final;

 private:
  static constexpr std::array<BinningValue, 3u> cylindricalBinning = {
      BinningValue::binZ, BinningValue::binR, BinningValue::binPhi};

  static constexpr std::array<BinningValue, 3u> cartesianBinning = {
      BinningValue::binX, BinningValue::binY, BinningValue::binZ};

  /// Configuration object
  Config m_cfg;

  /// Binning values from parsing
  std::array<BinningValue, 3u> m_binningValues = {};

  /// Private acces method to the logger
  const Logger& logger() const { return *m_logger; }

  /// logging instance
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Experimental
}  // namespace Acts
