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

  /// @brief Accessors for histograms (const reference)
  const Acts::Experimental::Histogram2& nRecoVsPt() const {
    return m_nRecoVsPt;
  }
  const Acts::Experimental::Histogram2& nTruthMatchedVsPt() const {
    return m_nTruthMatchedVsPt;
  }
  const Acts::Experimental::Histogram2& nFakeVsPt() const {
    return m_nFakeVsPt;
  }
  const Acts::Experimental::Histogram2& nRecoVsEta() const {
    return m_nRecoVsEta;
  }
  const Acts::Experimental::Histogram2& nTruthMatchedVsEta() const {
    return m_nTruthMatchedVsEta;
  }
  const Acts::Experimental::Histogram2& nFakeVsEta() const {
    return m_nFakeVsEta;
  }
  const Acts::Experimental::Efficiency1& fakeRatioVsPt() const {
    return m_fakeRatioVsPt;
  }
  const Acts::Experimental::Efficiency1& fakeRatioVsEta() const {
    return m_fakeRatioVsEta;
  }
  const Acts::Experimental::Efficiency1& fakeRatioVsPhi() const {
    return m_fakeRatioVsPhi;
  }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  Acts::Experimental::Histogram2 m_nRecoVsPt;
  Acts::Experimental::Histogram2 m_nTruthMatchedVsPt;
  Acts::Experimental::Histogram2 m_nFakeVsPt;
  Acts::Experimental::Histogram2 m_nRecoVsEta;
  Acts::Experimental::Histogram2 m_nTruthMatchedVsEta;
  Acts::Experimental::Histogram2 m_nFakeVsEta;
  Acts::Experimental::Efficiency1 m_fakeRatioVsPt;
  Acts::Experimental::Efficiency1 m_fakeRatioVsEta;
  Acts::Experimental::Efficiency1 m_fakeRatioVsPhi;
};

}  // namespace ActsExamples
