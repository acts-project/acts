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

/// Tools to make track quality plots.
class TrackQualityPlotTool {
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
    std::optional<Acts::Experimental::ProfileHistogram1> completeness_vs_pT;
    std::optional<Acts::Experimental::ProfileHistogram1> completeness_vs_eta;
    std::optional<Acts::Experimental::ProfileHistogram1> completeness_vs_phi;
    std::optional<Acts::Experimental::ProfileHistogram1> purity_vs_pT;
    std::optional<Acts::Experimental::ProfileHistogram1> purity_vs_eta;
    std::optional<Acts::Experimental::ProfileHistogram1> purity_vs_phi;
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

 private:
  /// The Config class
  Config m_cfg;
  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
