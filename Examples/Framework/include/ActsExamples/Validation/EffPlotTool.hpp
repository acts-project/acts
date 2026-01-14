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

  /// @brief Nested Cache struct
  struct Cache {
    /// Tracking efficiency vs eta
    Acts::Experimental::Efficiency1 trackEff_vs_eta;
    /// Tracking efficiency vs phi
    Acts::Experimental::Efficiency1 trackEff_vs_phi;
    /// Tracking efficiency vs pT
    Acts::Experimental::Efficiency1 trackEff_vs_pT;
    /// Tracking efficiency vs log pT
    Acts::Experimental::Efficiency1 trackEff_vs_LogPt;
    /// Tracking efficiency vs low pT
    Acts::Experimental::Efficiency1 trackEff_vs_LowPt;
    /// Tracking efficiency vs d0
    Acts::Experimental::Efficiency1 trackEff_vs_d0;
    /// Tracking efficiency vs z0
    Acts::Experimental::Efficiency1 trackEff_vs_z0;
    /// Tracking efficiency vs distance to the closest truth particle
    Acts::Experimental::Efficiency1 trackEff_vs_DeltaR;
    /// Tracking efficiency vs production radius
    Acts::Experimental::Efficiency1 trackEff_vs_prodR;

    /// Tracking efficiency vs eta and phi
    Acts::Experimental::Efficiency2 trackEff_vs_eta_phi;
    /// Tracking efficiency vs eta and pT
    Acts::Experimental::Efficiency2 trackEff_vs_eta_pt;

    /// Tracking efficiency vs eta in different pT ranges
    std::vector<Acts::Experimental::Efficiency1> trackEff_vs_eta_inPtRanges;
    /// Tracking efficiency vs pT in different abs(eta) ranges
    std::vector<Acts::Experimental::Efficiency1> trackEff_vs_pT_inAbsEtaRanges;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  EffPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the efficiency plots
  ///
  /// @param cache the cache for efficiency plots
  void book(Cache& cache) const;

  /// @brief fill efficiency plots
  ///
  /// @param gctx geometry context
  /// @param cache cache object for efficiency plots
  /// @param truthParticle the truth Particle
  /// @param deltaR the distance to the closest truth particle
  /// @param status the reconstruction status
  void fill(const Acts::GeometryContext& gctx, Cache& cache,
            const SimParticleState& truthParticle, double deltaR,
            bool status) const;

 private:
  /// The Config class
  Config m_cfg;
  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
