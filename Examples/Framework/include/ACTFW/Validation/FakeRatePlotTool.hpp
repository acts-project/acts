// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <memory>
#include <string>

#include "ACTFW/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsFatras/EventData/Particle.hpp"

namespace FW {

// Tools to make fake rate plots to show tracking fake rate.
class FakeRatePlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Eta", PlotHelpers::Binning("#eta", 40, -4, 4)},
        {"Phi", PlotHelpers::Binning("#phi", 100, -3.15, 3.15)},
        {"Pt", PlotHelpers::Binning("pT [GeV/c]", 20, 0, 100)},
        {"Num", PlotHelpers::Binning("N", 30, -0.5, 29.5)}};
  };

  /// @brief Nested Cache struct
  struct FakeRatePlotCache {
    TH1F* nRecoTracks;             ///< number of reco tracks
    TH1F* nTruthMatchedTracks;     ///< number of truth-matched reco tracks
    TH1F* nFakeTracks;             ///< number of fake tracks
    TEfficiency* fakerate_vs_eta;  ///< Tracking fake rate vs eta
    TEfficiency* fakerate_vs_phi;  ///< Tracking fake rate vs phi
    TEfficiency* fakerate_vs_pT;   ///< Tracking fake rate vs pT
    //@Todo: make duplication number plots with duplication plot tool
    TProfile* duplicationNum_vs_eta;  ///< Tracking duplication number vs eta
    TProfile* duplicationNum_vs_phi;  ///< Tracking duplication number vs phi
    TProfile* duplicationNum_vs_pT;   ///< Tracking duplication number vs pT
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  FakeRatePlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the fake rate plots
  /// @param fakeRatePlotCache the cache for fake rate plots
  void book(FakeRatePlotCache& fakeRatePlotCache) const;

  /// @brief fill fake rate plots
  ///
  /// @param fakeRatePlotCache cache object for fake rate plots
  /// @param truthParticle the truth Particle
  /// @param status the reconstruction status
  void fill(FakeRatePlotCache& fakeRatePlotCache,
            const ActsFatras::Particle& truthParticle, bool status) const;

  /// @brief fill number of reco/truth-matched/fake tracks for a single
  /// Multi-trajectory
  ///
  /// @param fakeRatePlotCache cache object for fake rate plots
  /// @param truthParticle the truth Particle
  /// @param nTruthMatchedTracks the number of truth-Matched tracks
  /// @param nFakeTracks the number of fake tracks
  void fill(FakeRatePlotCache& fakeRatePlotCache,
            const ActsFatras::Particle& truthParticle,
            size_t nTruthMatchedTracks, size_t nFakeTracks) const;

  /// @brief write the fake rate plots to file
  /// @param fakeRatePlotCache cache object for fake rate plots
  void write(const FakeRatePlotCache& fakeRatePlotCache) const;

  /// @brief delete the fake rate plots
  /// @param fakeRatePlotCache cache object for fake rate plots
  void clear(FakeRatePlotCache& fakeRatePlotCache) const;

 private:
  Config m_cfg;                                  ///< The Config class
  std::unique_ptr<const Acts::Logger> m_logger;  ///< The logging instance

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace FW
