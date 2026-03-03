// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <string>

namespace ActsExamples {

/// Tools to make track info plots to show tracking track info.
class TrackSummaryPlotTool {
 public:
  using AxisVariant = Acts::Experimental::AxisVariant;
  using BoostRegularAxis = Acts::Experimental::BoostRegularAxis;
  using ProfileHistogram1 = Acts::Experimental::ProfileHistogram1;

  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, AxisVariant> varBinning = {
        {"Eta", BoostRegularAxis(40, -4, 4, "#eta")},
        {"Phi", BoostRegularAxis(100, -3.15, 3.15, "#phi")},
        {"Pt", BoostRegularAxis(40, 0, 100, "pT [GeV/c]")},
        {"Num", BoostRegularAxis(30, -0.5, 29.5, "N")}};
    /// Optional prefix for histogram names
    std::string prefix;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  TrackSummaryPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief fill reco track info w.r.t. fitted track parameters
  ///
  /// @param fittedParameters fitted track parameters of this track
  /// @param nStates number of track states
  /// @param nMeasurements number of measurements
  /// @param nOutliers number of outliers
  /// @param nHoles number of holes
  /// @param nSharedHits number of shared hits
  void fill(const Acts::BoundTrackParameters& fittedParameters,
            std::size_t nStates, std::size_t nMeasurements,
            std::size_t nOutliers, std::size_t nHoles, std::size_t nSharedHits);

  /// @brief Accessor for profile histograms map (const reference)
  const std::map<std::string, ProfileHistogram1>& profiles() const {
    return m_profiles;
  }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  std::map<std::string, ProfileHistogram1> m_profiles;
};

}  // namespace ActsExamples
