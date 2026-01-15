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

/// Tools to make track quality plots.
class TrackQualityPlotTool {
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
  TrackQualityPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief fill track quality w.r.t. fitted track parameters
  ///
  /// @param fittedParameters fitted track parameters of this track
  /// @param completeness completeness of the track
  /// @param purity purity of the track
  void fill(const Acts::BoundTrackParameters& fittedParameters,
            double completeness, double purity);

  /// @brief Accessors for histograms (const reference)
  const Acts::Experimental::ProfileHistogram1& completenessVsPt() const {
    return m_completenessVsPt;
  }
  const Acts::Experimental::ProfileHistogram1& completenessVsEta() const {
    return m_completenessVsEta;
  }
  const Acts::Experimental::ProfileHistogram1& completenessVsPhi() const {
    return m_completenessVsPhi;
  }
  const Acts::Experimental::ProfileHistogram1& purityVsPt() const {
    return m_purityVsPt;
  }
  const Acts::Experimental::ProfileHistogram1& purityVsEta() const {
    return m_purityVsEta;
  }
  const Acts::Experimental::ProfileHistogram1& purityVsPhi() const {
    return m_purityVsPhi;
  }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  Acts::Experimental::ProfileHistogram1 m_completenessVsPt;
  Acts::Experimental::ProfileHistogram1 m_completenessVsEta;
  Acts::Experimental::ProfileHistogram1 m_completenessVsPhi;
  Acts::Experimental::ProfileHistogram1 m_purityVsPt;
  Acts::Experimental::ProfileHistogram1 m_purityVsEta;
  Acts::Experimental::ProfileHistogram1 m_purityVsPhi;
};

}  // namespace ActsExamples
