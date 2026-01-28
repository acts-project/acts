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
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"

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
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Eta", PlotHelpers::Binning::Uniform("#eta", 40, -3.0, 3.0)},
        {"Phi", PlotHelpers::Binning::Uniform("#phi", 100, -std::numbers::pi,
                                              std::numbers::pi)},
        {"Pt", PlotHelpers::Binning::Uniform("pT [GeV/c]", 40, 0, 100)},
        {"LogPt",
         PlotHelpers::Binning::Logarithmic("pT [GeV/c]", 11, 0.1, 100)},
        {"LowPt", PlotHelpers::Binning::Uniform("pT [GeV/c]", 40, 0, 2)},
        {"D0", PlotHelpers::Binning::Uniform("d_0 [mm]", 200, -200, 200)},
        {"Z0", PlotHelpers::Binning::Uniform("z_0 [mm]", 50, -200, 200)},
        {"DeltaR", PlotHelpers::Binning::Uniform("#Delta R", 100, 0, 0.3)},
        {"prodR", PlotHelpers::Binning::Uniform("prod_R [mm]", 100, 0, 200)}};

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
  EffPlotTool(const EffPlotTool&) = delete;
  EffPlotTool& operator=(const EffPlotTool&) = delete;
  EffPlotTool(EffPlotTool&&) noexcept = default;
  EffPlotTool& operator=(EffPlotTool&&) noexcept = default;
  ~EffPlotTool();

  /// @brief book the efficiency plots
  ///
  void book();

  /// @brief fill efficiency plots
  ///
  /// @param gctx geometry context
  /// @param truthParticle the truth Particle
  /// @param deltaR the distance to the closest truth particle
  /// @param status the reconstruction status
  void fill(const Acts::GeometryContext& gctx,
            const SimParticleState& truthParticle, double deltaR, bool status);

  /// @brief write the efficiency plots to file
  ///
  void write();

 private:
  struct Impl;

  /// The Config class
  Config m_cfg;
  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
  std::unique_ptr<Impl> m_impl;

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
