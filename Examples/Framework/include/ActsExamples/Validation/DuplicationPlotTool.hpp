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

using Acts::Experimental::AxisVariant;
using Acts::Experimental::BoostRegularAxis;

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
  DuplicationPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief fill duplication ratio w.r.t. fitted track parameters
  ///
  /// @param fittedParameters fitted track parameters of this track
  /// @param status the (truth-matched) reconstructed track is duplicated or not
  void fill(const Acts::BoundTrackParameters& fittedParameters, bool status);

  /// @brief fill number of duplicated tracks for a truth particle seed
  ///
  /// @param truthParticle the truth Particle
  /// @param nDuplicatedTracks the number of matched tracks
  void fill(const SimParticleState& truthParticle, std::size_t nMatchedTracks);

  /// @brief Accessors for histogram maps (const reference)
  const std::map<std::string, Acts::Experimental::ProfileHistogram1>& profiles()
      const {
    return m_profiles;
  }
  const std::map<std::string, Acts::Experimental::Efficiency1>& efficiencies()
      const {
    return m_efficiencies;
  }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  std::map<std::string, Acts::Experimental::ProfileHistogram1> m_profiles;
  std::map<std::string, Acts::Experimental::Efficiency1> m_efficiencies;
};

}  // namespace ActsExamples
