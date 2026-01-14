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
#include <optional>
#include <string>

namespace ActsExamples {

/// Tools to make track info plots to show tracking track info.
class TrackSummaryPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, Acts::Experimental::AxisVariant> varBinning = {
        {"Eta", Acts::Experimental::BoostRegularAxis(40, -4, 4, "#eta")},
        {"Phi", Acts::Experimental::BoostRegularAxis(100, -3.15, 3.15, "#phi")},
        {"Pt", Acts::Experimental::BoostRegularAxis(40, 0, 100, "pT [GeV/c]")},
        {"Num", Acts::Experimental::BoostRegularAxis(30, -0.5, 29.5, "N")}};
  };

  /// @brief Nested Cache struct
  struct Cache {
    /// Number of total states vs eta
    std::optional<Acts::Experimental::ProfileHistogram1> nStates_vs_eta;
    /// Number of non-outlier measurements vs eta
    std::optional<Acts::Experimental::ProfileHistogram1> nMeasurements_vs_eta;
    /// Number of holes vs eta
    std::optional<Acts::Experimental::ProfileHistogram1> nHoles_vs_eta;
    /// Number of outliers vs eta
    std::optional<Acts::Experimental::ProfileHistogram1> nOutliers_vs_eta;
    /// Number of Shared Hits vs eta
    std::optional<Acts::Experimental::ProfileHistogram1> nSharedHits_vs_eta;
    /// Number of total states vs pt
    std::optional<Acts::Experimental::ProfileHistogram1> nStates_vs_pt;
    /// Number of non-outlier measurements vs pt
    std::optional<Acts::Experimental::ProfileHistogram1> nMeasurements_vs_pt;
    /// Number of holes vs pt
    std::optional<Acts::Experimental::ProfileHistogram1> nHoles_vs_pt;
    /// Number of outliers vs pt
    std::optional<Acts::Experimental::ProfileHistogram1> nOutliers_vs_pt;
    /// Number of Shared Hits vs pt
    std::optional<Acts::Experimental::ProfileHistogram1> nSharedHits_vs_pt;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  TrackSummaryPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the track info plots
  ///
  /// @param cache the cache for track info plots
  /// @param prefix a prefix prepended to the name, concatenation with '_'
  void book(Cache& cache, const std::string& prefix = "") const;

  /// @brief fill reco track info w.r.t. fitted track parameters
  ///
  /// @param cache cache object for track info plots
  /// @param fittedParameters fitted track parameters of this track
  /// @param nStates number of track states
  /// @param nMeasurements number of measurements
  /// @param nOutliers number of outliers
  /// @param nHoles number of holes
  void fill(Cache& cache, const Acts::BoundTrackParameters& fittedParameters,
            std::size_t nStates, std::size_t nMeasurements,
            std::size_t Outliers, std::size_t nHoles,
            std::size_t nSharedHits) const;

 private:
  /// The Config class
  Config m_cfg;
  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
