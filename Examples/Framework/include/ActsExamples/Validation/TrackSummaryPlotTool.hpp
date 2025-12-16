// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"
#include "ActsExamples/Validation/BoostHistogramWrappers.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <string>

namespace ActsExamples {

/// Tools to make track info plots to show tracking track info.
class TrackSummaryPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Eta", PlotHelpers::Binning::Uniform("#eta", 40, -4, 4)},
        {"Phi", PlotHelpers::Binning::Uniform("#phi", 100, -3.15, 3.15)},
        {"Pt", PlotHelpers::Binning::Uniform("pT [GeV/c]", 40, 0, 100)},
        {"Num", PlotHelpers::Binning::Uniform("N", 30, -0.5, 29.5)}};
  };

  /// @brief Nested Cache struct
  struct Cache {
    /// Number of total states vs eta
    BoostProfileHistogram nStates_vs_eta;
    /// Number of non-outlier measurements vs eta
    BoostProfileHistogram nMeasurements_vs_eta;
    /// Number of holes vs eta
    BoostProfileHistogram nHoles_vs_eta;
    /// Number of outliers vs eta
    BoostProfileHistogram nOutliers_vs_eta;
    /// Number of Shared Hits vs eta
    BoostProfileHistogram nSharedHits_vs_eta;
    /// Number of total states vs pt
    BoostProfileHistogram nStates_vs_pt;
    /// Number of non-outlier measurements vs pt
    BoostProfileHistogram nMeasurements_vs_pt;
    /// Number of holes vs pt
    BoostProfileHistogram nHoles_vs_pt;
    /// Number of outliers vs pt
    BoostProfileHistogram nOutliers_vs_pt;
    /// Number of Shared Hits vs pt
    BoostProfileHistogram nSharedHits_vs_pt;
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
