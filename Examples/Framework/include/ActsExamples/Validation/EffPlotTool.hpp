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
  using AxisVariant = Acts::Experimental::AxisVariant;
  using BoostLogAxis = Acts::Experimental::BoostLogAxis;
  using BoostRegularAxis = Acts::Experimental::BoostRegularAxis;
  using Efficiency1 = Acts::Experimental::Efficiency1;
  using Efficiency2 = Acts::Experimental::Efficiency2;

  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, AxisVariant> varBinning = {
        {"Eta", BoostRegularAxis(40, -3.0, 3.0, "#eta")},
        {"Phi",
         BoostRegularAxis(100, -std::numbers::pi, std::numbers::pi, "#phi")},
        {"Pt", BoostRegularAxis(40, 0, 100, "pT [GeV/c]")},
        {"LogPt", BoostLogAxis(11, 0.1, 100, "pT [GeV/c]")},
        {"LowPt", BoostRegularAxis(40, 0, 2, "pT [GeV/c]")},
        {"D0", BoostRegularAxis(200, -200, 200, "d_0 [mm]")},
        {"Z0", BoostRegularAxis(50, -200, 200, "z_0 [mm]")},
        {"DeltaR", BoostRegularAxis(100, 0, 0.3, "#Delta R")},
        {"prodR", BoostRegularAxis(100, 0, 200, "prod_R [mm]")}};

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

  /// @brief Accessors for efficiency maps (const reference)
  const std::map<std::string, Efficiency1>& efficiencies1D() const {
    return m_efficiencies1D;
  }
  const std::map<std::string, Efficiency2>& efficiencies2D() const {
    return m_efficiencies2D;
  }
  const std::vector<Efficiency1>& trackEffVsEtaInPtRanges() const {
    return m_trackEffVsEtaInPtRanges;
  }
  const std::vector<Efficiency1>& trackEffVsPtInAbsEtaRanges() const {
    return m_trackEffVsPtInAbsEtaRanges;
  }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  std::map<std::string, Efficiency1> m_efficiencies1D;
  std::map<std::string, Efficiency2> m_efficiencies2D;
  std::vector<Efficiency1> m_trackEffVsEtaInPtRanges;
  std::vector<Efficiency1> m_trackEffVsPtInAbsEtaRanges;
};

}  // namespace ActsExamples
