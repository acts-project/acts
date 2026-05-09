// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <string>

namespace ActsExamples {

/// Tools to make track quality plots.
class TrackQualityPlotTool {
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
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  TrackQualityPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief fill track quality w.r.t. fitted track parameters
  ///
  /// @param fittedParameters fitted track parameters of this track
  /// @param completeness completeness of the track
  /// @param purity purity of the track
  void fill(const Acts::BoundTrackParameters& fittedParameters,
            double completeness, double purity);

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
