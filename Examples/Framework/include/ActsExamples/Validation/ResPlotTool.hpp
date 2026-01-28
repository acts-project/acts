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
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"

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
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Eta", PlotHelpers::Binning::Uniform("#eta", 40, -4, 4)},
        {"Pt", PlotHelpers::Binning::Uniform("pT [GeV/c]", 40, 0, 100)},
        {"Pull", PlotHelpers::Binning::Uniform("pull", 100, -5, 5)},
        {"Residual_d0",
         PlotHelpers::Binning::Uniform("r_{d0} [mm]", 100, -0.5, 0.5)},
        {"Residual_z0",
         PlotHelpers::Binning::Uniform("r_{z0} [mm]", 100, -0.5, 0.5)},
        {"Residual_phi",
         PlotHelpers::Binning::Uniform("r_{#phi} [rad]", 100, -0.01, 0.01)},
        {"Residual_theta",
         PlotHelpers::Binning::Uniform("r_{#theta} [rad]", 100, -0.01, 0.01)},
        {"Residual_qop",
         PlotHelpers::Binning::Uniform("r_{q/p} [c/GeV]", 100, -0.1, 0.1)},
        {"Residual_t",
         PlotHelpers::Binning::Uniform("r_{t} [s]", 100, -1000, 1000)}};
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  ResPlotTool(const Config& cfg, Acts::Logging::Level lvl);
  ResPlotTool(const ResPlotTool&) = delete;
  ResPlotTool& operator=(const ResPlotTool&) = delete;
  ResPlotTool(ResPlotTool&&) noexcept = default;
  ResPlotTool& operator=(ResPlotTool&&) noexcept = default;
  ~ResPlotTool();

  /// @brief book the histograms
  ///
  void book();

  /// @brief fill the histograms
  ///
  /// @param gctx the geometry context
  /// @param truthParticle the truth particle
  /// @param fittedParamters the fitted parameters at perigee surface
  void fill(const Acts::GeometryContext& gctx,
            const SimParticleState& truthParticle,
            const Acts::BoundTrackParameters& fittedParamters);

  /// @brief extract the details of the residual/pull plots and fill details
  ///
  /// into separate histograms
  void refinement();

  /// @brief write the histograms to output file
  ///
  void write();

 private:
  struct Impl;

  /// The config class
  Config m_cfg;
  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
  std::unique_ptr<Impl> m_impl;

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
