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

#include <cstddef>
#include <map>
#include <memory>
#include <string>

namespace ActsExamples {

using Acts::Experimental::AxisVariant;

/// Tools to make track info plots to show tracking track info.
class TrackSummaryPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, AxisVariant> varBinning = {
        {"Eta", Acts::Experimental::BoostRegularAxis(40, -4, 4, "#eta")},
        {"Phi", Acts::Experimental::BoostRegularAxis(100, -3.15, 3.15, "#phi")},
        {"Pt", Acts::Experimental::BoostRegularAxis(40, 0, 100, "pT [GeV/c]")},
        {"Num", Acts::Experimental::BoostRegularAxis(30, -0.5, 29.5, "N")}};
    /// Optional prefix for histogram names
    std::string prefix;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  TrackSummaryPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief fill reco track info w.r.t. fitted track parameters
  ///
  /// @param fittedParameters fitted track parameters of this track
  /// @param nStates number of track states
  /// @param nMeasurements number of measurements
  /// @param nOutliers number of outliers
  /// @param nHoles number of holes
  /// @param nSharedHits number of shared hits
  void fill(const Acts::BoundTrackParameters& fittedParameters,
            std::size_t nStates, std::size_t nMeasurements,
            std::size_t nOutliers, std::size_t nHoles, std::size_t nSharedHits);

  /// @brief Accessors for histograms (const reference)
  const Acts::Experimental::ProfileHistogram1& nStatesVsEta() const {
    return m_nStatesVsEta;
  }
  const Acts::Experimental::ProfileHistogram1& nMeasurementsVsEta() const {
    return m_nMeasurementsVsEta;
  }
  const Acts::Experimental::ProfileHistogram1& nHolesVsEta() const {
    return m_nHolesVsEta;
  }
  const Acts::Experimental::ProfileHistogram1& nOutliersVsEta() const {
    return m_nOutliersVsEta;
  }
  const Acts::Experimental::ProfileHistogram1& nSharedHitsVsEta() const {
    return m_nSharedHitsVsEta;
  }
  const Acts::Experimental::ProfileHistogram1& nStatesVsPt() const {
    return m_nStatesVsPt;
  }
  const Acts::Experimental::ProfileHistogram1& nMeasurementsVsPt() const {
    return m_nMeasurementsVsPt;
  }
  const Acts::Experimental::ProfileHistogram1& nHolesVsPt() const {
    return m_nHolesVsPt;
  }
  const Acts::Experimental::ProfileHistogram1& nOutliersVsPt() const {
    return m_nOutliersVsPt;
  }
  const Acts::Experimental::ProfileHistogram1& nSharedHitsVsPt() const {
    return m_nSharedHitsVsPt;
  }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  Acts::Experimental::ProfileHistogram1 m_nStatesVsEta;
  Acts::Experimental::ProfileHistogram1 m_nMeasurementsVsEta;
  Acts::Experimental::ProfileHistogram1 m_nHolesVsEta;
  Acts::Experimental::ProfileHistogram1 m_nOutliersVsEta;
  Acts::Experimental::ProfileHistogram1 m_nSharedHitsVsEta;
  Acts::Experimental::ProfileHistogram1 m_nStatesVsPt;
  Acts::Experimental::ProfileHistogram1 m_nMeasurementsVsPt;
  Acts::Experimental::ProfileHistogram1 m_nHolesVsPt;
  Acts::Experimental::ProfileHistogram1 m_nOutliersVsPt;
  Acts::Experimental::ProfileHistogram1 m_nSharedHitsVsPt;
};

}  // namespace ActsExamples
