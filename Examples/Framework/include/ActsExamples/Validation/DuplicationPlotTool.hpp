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

/// Tools to make duplication ratio (formerly called duplication rate) and
/// duplication number plots to show tracking duplication.
///
/// The duplication ratio is evaluated on truth-matched reco tracks. If there
/// is more than one track matched to the same truth particle, the reco track
/// with the highest matching probability is tagged as 'real' and the others are
/// 'duplicated'.
class DuplicationPlotTool {
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
    /// Number of duplicated tracks vs pT
    Acts::Experimental::ProfileHistogram1 nDuplicated_vs_pT;
    /// Number of duplicated tracks vs eta
    Acts::Experimental::ProfileHistogram1 nDuplicated_vs_eta;
    /// Number of duplicated tracks vs phi
    Acts::Experimental::ProfileHistogram1 nDuplicated_vs_phi;
    /// Tracking duplication ratio vs pT
    Acts::Experimental::Efficiency1 duplicationRatio_vs_pT;
    /// Tracking duplication ratio vs eta
    Acts::Experimental::Efficiency1 duplicationRatio_vs_eta;
    /// Tracking duplication ratio vs phi
    Acts::Experimental::Efficiency1 duplicationRatio_vs_phi;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  DuplicationPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the duplication plots
  ///
  /// @param cache the cache for duplication plots
  void book(Cache& cache) const;

  /// @brief fill duplication ratio w.r.t. fitted track parameters
  ///
  /// @param cache cache object for duplication plots
  /// @param fittedParameters fitted track parameters of this track
  /// @param status the (truth-matched) reconstructed track is duplicated or not
  void fill(Cache& cache, const Acts::BoundTrackParameters& fittedParameters,
            bool status) const;

  /// @brief fill number of duplicated tracks for a truth particle seed
  ///
  /// @param cache cache object for duplication plots
  /// @param truthParticle the truth Particle
  /// @param nDuplicatedTracks the number of duplicated tracks
  void fill(Cache& cache, const SimParticleState& truthParticle,
            std::size_t nDuplicatedTracks) const;

 private:
  /// The Config class
  Config m_cfg;
  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
