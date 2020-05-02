// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <memory>
#include <string>

#include "ACTFW/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsFatras/EventData/Particle.hpp"

namespace FW {

// Tools to make duplication number&rate plots to show tracking duplication.
class DuplicationPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Eta", PlotHelpers::Binning("#eta", 40, -4, 4)},
        {"Phi", PlotHelpers::Binning("#phi", 100, -3.15, 3.15)},
        {"Pt", PlotHelpers::Binning("pT [GeV/c]", 20, 0, 100)},
        {"Num", PlotHelpers::Binning("N", 30, -0.5, 29.5)}};
  };

  /// @brief Nested Cache struct
  struct DuplicationPlotCache {
    TProfile* duplicationNum_vs_eta;  ///< Tracking duplication number vs eta
    TProfile* duplicationNum_vs_phi;  ///< Tracking duplication number vs phi
    TProfile* duplicationNum_vs_pT;   ///< Tracking duplication number vs pT
    TEfficiency* duplicationRate_vs_eta;  ///< Tracking duplication rate vs eta
    TEfficiency* duplicationRate_vs_phi;  ///< Tracking duplication rate vs phi
    TEfficiency* duplicationRate_vs_pT;   ///< Tracking duplication rate vs pT
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  DuplicationPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the duplication plots
  /// @param duplicationPlotCache the cache for duplication plots
  void book(DuplicationPlotCache& duplicationPlotCache) const;

  /// @brief fill duplication plots
  ///
  /// @param duplicationPlotCache cache object for duplication plots
  /// @param truthParticle the truth Particle
  /// @param status the reconstruction status
  void fill(DuplicationPlotCache& duplicationPlotCache,
            const ActsFatras::Particle& truthParticle, bool status) const;

  /// @brief write the duplication plots to file
  /// @param duplicationPlotCache cache object for duplication plots
  void write(const DuplicationPlotCache& duplicationPlotCache) const;

  /// @brief delete the duplication plots
  /// @param duplicationPlotCache cache object for duplication plots
  void clear(DuplicationPlotCache& duplicationPlotCache) const;

 private:
  Config m_cfg;                                  ///< The Config class
  std::unique_ptr<const Acts::Logger> m_logger;  ///< The logging instance

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace FW
