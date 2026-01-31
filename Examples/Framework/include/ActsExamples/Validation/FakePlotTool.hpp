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
#include "ActsExamples/EventData/SimParticle.hpp"

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
  using AxisVariant = Acts::Experimental::AxisVariant;
  using BoostRegularAxis = Acts::Experimental::BoostRegularAxis;
  using Efficiency1 = Acts::Experimental::Efficiency1;
  using Histogram2 = Acts::Experimental::Histogram2;

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
  FakePlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief fill fake ratio w.r.t. fitted track parameters
  ///
  /// @param fittedParameters fitted track parameters of this track
  /// @param status the reconstructed track is fake or not
  void fill(const Acts::BoundTrackParameters& fittedParameters, bool status);

  /// @brief fill number of reco/truth-matched/fake tracks for a truth particle
  ///
  /// @param truthParticle the truth Particle
  /// @param nTruthMatchedTracks the number of truth-matched tracks
  /// @param nFakeTracks the number of fake tracks
  void fill(const SimParticleState& truthParticle,
            std::size_t nTruthMatchedTracks, std::size_t nFakeTracks);

  /// @brief Accessors for histogram maps (const reference)
  const std::map<std::string, Histogram2>& histograms() const {
    return m_histograms;
  }
  const std::map<std::string, Efficiency1>& efficiencies() const {
    return m_efficiencies;
  }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  std::map<std::string, Histogram2> m_histograms;
  std::map<std::string, Efficiency1> m_efficiencies;
};

}  // namespace ActsExamples
