// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace ActsExamples {

/// Tools to make efficiency plots to show tracking efficiency.
/// For the moment, the efficiency is taken as the fraction of successfully
/// smoothed track over all tracks
class EffPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, Acts::Experimental::AxisVariant> varBinning = {
        {"Eta", Acts::Experimental::BoostRegularAxis(40, -3.0, 3.0, "#eta")},
        {"Phi", Acts::Experimental::BoostRegularAxis(100, -std::numbers::pi,
                                                      std::numbers::pi, "#phi")},
        {"Pt", Acts::Experimental::BoostRegularAxis(40, 0, 100, "pT [GeV/c]")},
        {"LogPt", Acts::Experimental::BoostLogAxis(11, 0.1, 100, "pT [GeV/c]")},
        {"LowPt", Acts::Experimental::BoostRegularAxis(40, 0, 2, "pT [GeV/c]")},
        {"D0", Acts::Experimental::BoostRegularAxis(200, -200, 200, "d_0 [mm]")},
        {"Z0", Acts::Experimental::BoostRegularAxis(50, -200, 200, "z_0 [mm]")},
        {"DeltaR", Acts::Experimental::BoostRegularAxis(100, 0, 0.3, "#Delta R")},
        {"prodR",
         Acts::Experimental::BoostRegularAxis(100, 0, 200, "prod_R [mm]")}};

    double minTruthPt = 1.0 * Acts::UnitConstants::GeV;

    std::vector<std::pair<double, double>> truthPtRangesForEta = {
        {0.9 * Acts::UnitConstants::GeV, 1.1 * Acts::UnitConstants::GeV},
        {9 * Acts::UnitConstants::GeV, 11 * Acts::UnitConstants::GeV},
        {90 * Acts::UnitConstants::GeV, 110 * Acts::UnitConstants::GeV}};

    std::vector<std::pair<double, double>> truthAbsEtaRangesForPt = {
        {0, 0.2}, {0, 0.8}, {1, 2}, {2, 3}};

    /// Beamline to estimate d0 and z0
    std::shared_ptr<Acts::Surface> beamline =
        Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  EffPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief fill efficiency plots
  ///
  /// @param gctx geometry context
  /// @param truthParticle the truth Particle
  /// @param deltaR the distance to the closest truth particle
  /// @param status the reconstruction status
  void fill(const Acts::GeometryContext& gctx,
            const SimParticleState& truthParticle, double deltaR, bool status);

  /// @brief Accessors for histograms (const reference)
  const Acts::Experimental::Efficiency1& trackEffVsEta() const {
    return m_trackEffVsEta;
  }
  const Acts::Experimental::Efficiency1& trackEffVsPhi() const {
    return m_trackEffVsPhi;
  }
  const Acts::Experimental::Efficiency1& trackEffVsPt() const {
    return m_trackEffVsPt;
  }
  const Acts::Experimental::Efficiency1& trackEffVsLogPt() const {
    return m_trackEffVsLogPt;
  }
  const Acts::Experimental::Efficiency1& trackEffVsLowPt() const {
    return m_trackEffVsLowPt;
  }
  const Acts::Experimental::Efficiency1& trackEffVsD0() const {
    return m_trackEffVsD0;
  }
  const Acts::Experimental::Efficiency1& trackEffVsZ0() const {
    return m_trackEffVsZ0;
  }
  const Acts::Experimental::Efficiency1& trackEffVsDeltaR() const {
    return m_trackEffVsDeltaR;
  }
  const Acts::Experimental::Efficiency1& trackEffVsProdR() const {
    return m_trackEffVsProdR;
  }
  const Acts::Experimental::Efficiency2& trackEffVsEtaPhi() const {
    return m_trackEffVsEtaPhi;
  }
  const Acts::Experimental::Efficiency2& trackEffVsEtaPt() const {
    return m_trackEffVsEtaPt;
  }
  const std::vector<Acts::Experimental::Efficiency1>& trackEffVsEtaInPtRanges()
      const {
    return m_trackEffVsEtaInPtRanges;
  }
  const std::vector<Acts::Experimental::Efficiency1>&
  trackEffVsPtInAbsEtaRanges() const {
    return m_trackEffVsPtInAbsEtaRanges;
  }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  Acts::Experimental::Efficiency1 m_trackEffVsEta;
  Acts::Experimental::Efficiency1 m_trackEffVsPhi;
  Acts::Experimental::Efficiency1 m_trackEffVsPt;
  Acts::Experimental::Efficiency1 m_trackEffVsLogPt;
  Acts::Experimental::Efficiency1 m_trackEffVsLowPt;
  Acts::Experimental::Efficiency1 m_trackEffVsD0;
  Acts::Experimental::Efficiency1 m_trackEffVsZ0;
  Acts::Experimental::Efficiency1 m_trackEffVsDeltaR;
  Acts::Experimental::Efficiency1 m_trackEffVsProdR;
  Acts::Experimental::Efficiency2 m_trackEffVsEtaPhi;
  Acts::Experimental::Efficiency2 m_trackEffVsEtaPt;
  std::vector<Acts::Experimental::Efficiency1> m_trackEffVsEtaInPtRanges;
  std::vector<Acts::Experimental::Efficiency1> m_trackEffVsPtInAbsEtaRanges;
};

}  // namespace ActsExamples
