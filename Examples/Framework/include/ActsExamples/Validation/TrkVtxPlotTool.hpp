// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <map>
#include <memory>
#include <string>

class TEfficiency;
namespace ActsFatras {
class Particle;
}  // namespace ActsFatras

namespace ActsExamples {

// Tools to make efficiency plots to show tracking efficiency.
// For the moment, the efficiency is taken as the fraction of successfully
// smoothed track over all tracks
class TrkVtxPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Zgen", PlotHelpers::Binning("z_{gen} (mm)", 70, -65, 5)},
        {"Zrec", PlotHelpers::Binning("z_{rec} (mm)", 70, -65, 5)},
        {"Tar", PlotHelpers::Binning("z_{gen} (mm)", 5, -60, 6)},
        {"Res", PlotHelpers::Binning("z_{rec}-z_{gen} (mm)", 1000, -50, 50)},
    };
  };

  /// @brief Nested Cache struct
  struct TrkVtxPlotCache {
    TH1F* resVtxz{nullptr};   ///< Tracking efficiency vs pT
    TEfficiency* eff_vs_zgen{nullptr};   ///< Tracking efficiency vs pT
    TH2F* resVtxz_vs_zgen{nullptr};  ///< Tracking efficiency vs eta
    TH2F* zrec_vs_zgen{nullptr};  ///< Tracking efficiency vs eta
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  TrkVtxPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the efficiency plots
  ///
  /// @param trkVtxPlotCache the cache for efficiency plots
  void book(TrkVtxPlotCache& trkVtxPlotCache) const;

  /// @brief fill efficiency plots
  ///
  /// @param trkVtxPlotCache cache object for efficiency plots
  /// @param truthParticle the truth Particle
  /// @param status the reconstruction status
  void fill(TrkVtxPlotTool::TrkVtxPlotCache& trkVtxPlotCache, double zrec, double zgen) const;

  /// @brief write the efficiency plots to file
  ///
  /// @param trkVtxPlotCache cache object for efficiency plots
  void write(const TrkVtxPlotCache& trkVtxPlotCache) const;

  /// @brief delete the efficiency plots
  ///
  /// @param trkVtxPlotCache cache object for efficiency plots
  void clear(TrkVtxPlotCache& trkVtxPlotCache) const;

 private:
  Config m_cfg;                                  ///< The Config class
  std::unique_ptr<const Acts::Logger> m_logger;  ///< The logging instance

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
