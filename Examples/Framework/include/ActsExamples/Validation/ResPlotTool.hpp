// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace ActsExamples {

/// Tools to make hists to show residual, i.e. smoothed_parameter -
/// truth_parameter, and pull, i.e. (smoothed_parameter -
/// truth_parameter)/smoothed_paramter_error, of track parameters at perigee
/// surface
class ResPlotTool {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// parameter sets to do plots
    std::vector<std::string> paramNames = {"d0",    "z0",  "phi",
                                           "theta", "qop", "t"};

    /// Binning info for variables
    std::map<std::string, Acts::Experimental::AxisVariant> varBinning = {
        {"Eta", Acts::Experimental::BoostRegularAxis(40, -4, 4, "#eta")},
        {"Pt", Acts::Experimental::BoostRegularAxis(40, 0, 100, "pT [GeV/c]")},
        {"Pull", Acts::Experimental::BoostRegularAxis(100, -5, 5, "pull")},
        {"Residual_d0",
         Acts::Experimental::BoostRegularAxis(100, -0.5, 0.5, "r_{d0} [mm]")},
        {"Residual_z0",
         Acts::Experimental::BoostRegularAxis(100, -0.5, 0.5, "r_{z0} [mm]")},
        {"Residual_phi", Acts::Experimental::BoostRegularAxis(
                             100, -0.01, 0.01, "r_{#phi} [rad]")},
        {"Residual_theta", Acts::Experimental::BoostRegularAxis(
                               100, -0.01, 0.01, "r_{#theta} [rad]")},
        {"Residual_qop", Acts::Experimental::BoostRegularAxis(
                             100, -0.1, 0.1, "r_{q/p} [c/GeV]")},
        {"Residual_t",
         Acts::Experimental::BoostRegularAxis(100, -1000, 1000, "r_{t} [s]")}};
  };

  /// @brief Nested Cache struct
  struct Cache {
    /// Residual distribution
    std::map<std::string, Acts::Experimental::Histogram1> res;
    /// Residual vs eta scatter plot
    std::map<std::string, Acts::Experimental::Histogram2> res_vs_eta;
    /// Residual vs pT scatter plot
    std::map<std::string, Acts::Experimental::Histogram2> res_vs_pT;

    /// Pull distribution
    std::map<std::string, Acts::Experimental::Histogram1> pull;
    /// Pull vs eta scatter plot
    std::map<std::string, Acts::Experimental::Histogram2> pull_vs_eta;
    /// Pull vs pT scatter plot
    std::map<std::string, Acts::Experimental::Histogram2> pull_vs_pT;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  ResPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the histograms
  ///
  /// @param cache the cache for residual/pull histograms
  void book(Cache& cache) const;

  /// @brief fill the histograms
  ///
  /// @param cache the cache for residual/pull histograms
  /// @param gctx the geometry context
  /// @param truthParticle the truth particle
  /// @param fittedParamters the fitted parameters at perigee surface
  void fill(Cache& cache, const Acts::GeometryContext& gctx,
            const SimParticleState& truthParticle,
            const Acts::BoundTrackParameters& fittedParamters) const;

 private:
  /// The config class
  Config m_cfg;
  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
