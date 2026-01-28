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
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <string>

namespace ActsExamples {

/// Tools to make fake ratio plots.
///
/// The fake ratio (formerly called fake rate) is evaluated on all reco tracks.
/// A track is 'fake' if it's not matched with truth.
class FakePlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Eta", PlotHelpers::Binning::Uniform("#eta", 40, -4, 4)},
        {"Phi", PlotHelpers::Binning::Uniform("#phi", 100, -3.15, 3.15)},
        {"Pt", PlotHelpers::Binning::Uniform("pT [GeV/c]", 40, 0, 100)},
        {"Num", PlotHelpers::Binning::Uniform("N", 30, -0.5, 29.5)}};
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  FakePlotTool(const Config& cfg, Acts::Logging::Level lvl);
  FakePlotTool(const FakePlotTool&) = delete;
  FakePlotTool& operator=(const FakePlotTool&) = delete;
  FakePlotTool(FakePlotTool&&) noexcept = default;
  FakePlotTool& operator=(FakePlotTool&&) noexcept = default;
  ~FakePlotTool();

  /// @brief book the fake ratio plots
  ///
  void book();

  /// @brief fill fake ratio w.r.t. fitted track parameters
  ///
  /// @param fittedParameters fitted track parameters of this track
  /// @param status the reconstructed track is fake or not
  void fill(const Acts::BoundTrackParameters& fittedParameters, bool status);

  /// @brief fill number of reco/truth-matched/fake tracks for a truth particle
  /// seed
  ///
  /// @param truthParticle the truth Particle
  /// @param nTruthMatchedTracks the number of truth-matched tracks
  /// @param nFakeTracks the number of fake tracks
  void fill(const SimParticleState& truthParticle,
            std::size_t nTruthMatchedTracks, std::size_t nFakeTracks);

  /// @brief write the fake ratio plots to file
  ///
  void write();

 private:
  struct Impl;

  /// The Config class
  Config m_cfg;
  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
  std::unique_ptr<Impl> m_impl;

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
