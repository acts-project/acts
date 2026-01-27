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
  void fill(const SimParticleState& truthParticle,
            std::size_t nMatchedTracks);

  /// @brief Accessors for histograms (const reference)
  const Acts::Experimental::ProfileHistogram1& nDuplicatedVsPt() const {
    return m_nDuplicatedVsPt;
  }
  const Acts::Experimental::ProfileHistogram1& nDuplicatedVsEta() const {
    return m_nDuplicatedVsEta;
  }
  const Acts::Experimental::ProfileHistogram1& nDuplicatedVsPhi() const {
    return m_nDuplicatedVsPhi;
  }
  const Acts::Experimental::Efficiency1& duplicationRatioVsPt() const {
    return m_duplicationRatioVsPt;
  }
  const Acts::Experimental::Efficiency1& duplicationRatioVsEta() const {
    return m_duplicationRatioVsEta;
  }
  const Acts::Experimental::Efficiency1& duplicationRatioVsPhi() const {
    return m_duplicationRatioVsPhi;
  }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  Acts::Experimental::ProfileHistogram1 m_nDuplicatedVsPt;
  Acts::Experimental::ProfileHistogram1 m_nDuplicatedVsEta;
  Acts::Experimental::ProfileHistogram1 m_nDuplicatedVsPhi;
  Acts::Experimental::Efficiency1 m_duplicationRatioVsPt;
  Acts::Experimental::Efficiency1 m_duplicationRatioVsEta;
  Acts::Experimental::Efficiency1 m_duplicationRatioVsPhi;
};

}  // namespace ActsExamples
