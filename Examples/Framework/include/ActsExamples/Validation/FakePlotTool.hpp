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

class TEfficiency;
class TH2F;

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
        {"Eta", PlotHelpers::Binning("#eta", 40, -4, 4)},
        {"Phi", PlotHelpers::Binning("#phi", 100, -3.15, 3.15)},
        {"Pt", PlotHelpers::Binning("pT [GeV/c]", 40, 0, 100)},
        {"Num", PlotHelpers::Binning("N", 30, -0.5, 29.5)}};
  };

  /// @brief Nested Cache struct
  struct Cache {
    /// Number of reco tracks vs pT scatter plot
    TH2F* nReco_vs_pT;
    /// Number of truth-matched reco tracks vs pT scatter plot
    TH2F* nTruthMatched_vs_pT;
    /// Number of fake (truth-unmatched) tracks vs pT scatter plot
    TH2F* nFake_vs_pT;
    /// Number of reco tracks vs eta scatter plot
    TH2F* nReco_vs_eta;
    /// Number of truth-matched reco tracks vs eta scatter plot
    TH2F* nTruthMatched_vs_eta;
    /// Number of fake (truth-unmatched) tracks vs eta scatter plot
    TH2F* nFake_vs_eta;
    /// Tracking fake ratio vs pT
    TEfficiency* fakeRatio_vs_pT;
    /// Tracking fake ratio vs eta
    TEfficiency* fakeRatio_vs_eta;
    /// Tracking fake ratio vs phi
    TEfficiency* fakeRatio_vs_phi;
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

  /// @brief write the fake ratio plots to file
  ///
  /// @param cache cache object for fake ratio plots
  void write(const Cache& cache) const;

  /// @brief delete the fake ratio plots
  ///
  /// @param cache cache object for fake ratio plots
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
