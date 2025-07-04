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

#include <cstddef>
#include <map>
#include <memory>
#include <string>

class TProfile;

namespace ActsExamples {

/// Tools to make track quality plots.
class TrackQualityPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Eta", PlotHelpers::Binning("#eta", 40, -4, 4)},
        {"Phi", PlotHelpers::Binning("#phi", 100, -3.15, 3.15)},
        {"Pt", PlotHelpers::Binning("pT [GeV/c]", 40, 0, 100)},
        {"Num", PlotHelpers::Binning("N", 30, -0.5, 29.5)}};
  };

  /// @brief Nested Cache struct
  struct Cache {
    TProfile* completeness_vs_pT;
    TProfile* completeness_vs_eta;
    TProfile* completeness_vs_phi;
    TProfile* purity_vs_pT;
    TProfile* purity_vs_eta;
    TProfile* purity_vs_phi;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  TrackQualityPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the track quality plots
  ///
  /// @param cache the cache for track quality plots
  void book(Cache& cache) const;

  /// @brief fill track quality w.r.t. fitted track parameters
  ///
  /// @param cache cache object for track quality plots
  /// @param fittedParameters fitted track parameters of this track
  /// @param completeness completeness of the track
  /// @param purity purity of the track
  void fill(Cache& cache, const Acts::BoundTrackParameters& fittedParameters,
            double completeness, double purity) const;

  /// @brief write the track quality plots to file
  ///
  /// @param cache cache object for track quality plots
  void write(const Cache& cache) const;

  /// @brief delete the track quality plots
  ///
  /// @param cache cache object for track quality plots
  void clear(Cache& cache) const;

 private:
  /// The Config class
  Config m_cfg;
  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
