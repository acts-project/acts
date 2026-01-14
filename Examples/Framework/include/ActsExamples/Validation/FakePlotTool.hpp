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
#include <optional>
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
    std::map<std::string, Acts::Experimental::AxisVariant> varBinning = {
        {"Eta", Acts::Experimental::BoostRegularAxis(40, -4, 4, "#eta")},
        {"Phi", Acts::Experimental::BoostRegularAxis(100, -3.15, 3.15, "#phi")},
        {"Pt", Acts::Experimental::BoostRegularAxis(40, 0, 100, "pT [GeV/c]")},
        {"Num", Acts::Experimental::BoostRegularAxis(30, -0.5, 29.5, "N")}};
  };

  /// @brief Nested Cache struct
  struct Cache {
    /// Number of reco tracks vs pT scatter plot
    Acts::Experimental::Histogram2 nReco_vs_pT;
    /// Number of truth-matched reco tracks vs pT scatter plot
    Acts::Experimental::Histogram2 nTruthMatched_vs_pT;
    /// Number of fake (truth-unmatched) tracks vs pT scatter plot
    Acts::Experimental::Histogram2 nFake_vs_pT;
    /// Number of reco tracks vs eta scatter plot
    Acts::Experimental::Histogram2 nReco_vs_eta;
    /// Number of truth-matched reco tracks vs eta scatter plot
    Acts::Experimental::Histogram2 nTruthMatched_vs_eta;
    /// Number of fake (truth-unmatched) tracks vs eta scatter plot
    Acts::Experimental::Histogram2 nFake_vs_eta;
    /// Tracking fake ratio vs pT
    Acts::Experimental::Efficiency1 fakeRatio_vs_pT;
    /// Tracking fake ratio vs eta
    Acts::Experimental::Efficiency1 fakeRatio_vs_eta;
    /// Tracking fake ratio vs phi
    Acts::Experimental::Efficiency1 fakeRatio_vs_phi;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  FakePlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the fake ratio plots
  ///
  /// @param cache the cache for fake ratio plots
  void book(Cache& cache) const;

  /// @brief fill fake ratio w.r.t. fitted track parameters
  ///
  /// @param cache cache object for fake ratio plots
  /// @param fittedParameters fitted track parameters of this track
  /// @param status the reconstructed track is fake or not
  void fill(Cache& cache, const Acts::BoundTrackParameters& fittedParameters,
            bool status) const;

  /// @brief fill number of reco/truth-matched/fake tracks for a truth particle
  /// seed
  ///
  /// @param cache cache object for fake ratio plots
  /// @param truthParticle the truth Particle
  /// @param nTruthMatchedTracks the number of truth-matched tracks
  /// @param nFakeTracks the number of fake tracks
  void fill(Cache& cache, const SimParticleState& truthParticle,
            std::size_t nTruthMatchedTracks, std::size_t nFakeTracks) const;

 private:
  /// The Config class
  Config m_cfg;
  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
