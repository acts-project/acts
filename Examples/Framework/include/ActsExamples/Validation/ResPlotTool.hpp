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
  using AxisVariant = Acts::Experimental::AxisVariant;
  using BoostRegularAxis = Acts::Experimental::BoostRegularAxis;
  using Histogram1 = Acts::Experimental::Histogram1;
  using Histogram2 = Acts::Experimental::Histogram2;
  using Histogram3 = Acts::Experimental::Histogram3;

  /// @brief Nested configuration struct
  struct Config {
    /// Track parameter names
    std::vector<std::string> paramNames = {"d0",    "z0",  "phi",
                                           "theta", "qop", "t"};

    std::string qOverPtName = "qOverPt";
    std::string relQoverPtName = "relQoverPt";

    /// Binning info for variables
    std::map<std::string, AxisVariant> varBinning = {
        {"Eta", BoostRegularAxis(40, -4, 4, "#eta")},
        {"Phi", BoostRegularAxis(40, -std::numbers::pi, std::numbers::pi,
                                 "#phi [rad]")},
        {"Pt", BoostRegularAxis(40, 0, 100, "pT [GeV/c]")},
        {"Pull", BoostRegularAxis(100, -5, 5, "pull")},
        {"Residual_d0", BoostRegularAxis(100, -0.5, 0.5, "r_{d0} [mm]")},
        {"Residual_z0", BoostRegularAxis(100, -0.5, 0.5, "r_{z0} [mm]")},
        {"Residual_phi", BoostRegularAxis(100, -0.01, 0.01, "r_{#phi} [rad]")},
        {"Residual_theta",
         BoostRegularAxis(100, -0.01, 0.01, "r_{#theta} [rad]")},
        {"Residual_qop", BoostRegularAxis(100, -0.1, 0.1, "r_{q/p} [c/GeV]")},
        {"Residual_t", BoostRegularAxis(100, -100, 100, "r_{t} [mm/c]")},
        {"Residual_qOverPt",
         BoostRegularAxis(100, -0.1, 0.1, "r_{q/pT} [c/GeV]")},
        {"Residual_relQoverPt",
         BoostRegularAxis(100, -0.1, 0.1, "r_{rel q/pT} [%]")}};
  };

  /// @param cfg Configuration struct
  /// @param level Message level declaration
  ResPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @param gctx the geometry context
  /// @param truthParticle the truth particle
  /// @param fittedParamters the fitted parameters at perigee surface
  void fill(const Acts::GeometryContext& gctx,
            const SimParticleState& truthParticle,
            const Acts::BoundTrackParameters& fittedParamters);

  const std::map<std::string, Histogram1>& res() const { return m_res; }
  const std::map<std::string, Histogram2>& resVsEta() const {
    return m_resVsEta;
  }
  const std::map<std::string, Histogram2>& resVsPt() const { return m_resVsPt; }
  const std::map<std::string, Histogram3>& resVsEtaPhi() const {
    return m_resVsEtaPhi;
  }
  const std::map<std::string, Histogram3>& resVsEtaPt() const {
    return m_resVsEtaPt;
  }
  const std::map<std::string, Histogram1>& pull() const { return m_pull; }
  const std::map<std::string, Histogram2>& pullVsEta() const {
    return m_pullVsEta;
  }
  const std::map<std::string, Histogram2>& pullVsPt() const {
    return m_pullVsPt;
  }
  const std::map<std::string, Histogram3>& pullVsEtaPhi() const {
    return m_pullVsEtaPhi;
  }
  const std::map<std::string, Histogram3>& pullVsEtaPt() const {
    return m_pullVsEtaPt;
  }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;

  std::unique_ptr<const Acts::Logger> m_logger;

  /// Residual distribution
  std::map<std::string, Histogram1> m_res;
  /// Residual vs eta scatter plot
  std::map<std::string, Histogram2> m_resVsEta;
  /// Residual vs pT scatter plot
  std::map<std::string, Histogram2> m_resVsPt;
  /// Residual vs eta-phi scatter plot
  std::map<std::string, Histogram3> m_resVsEtaPhi;
  /// Residual vs eta-pT scatter plot
  std::map<std::string, Histogram3> m_resVsEtaPt;

  /// Pull distribution
  std::map<std::string, Histogram1> m_pull;
  /// Pull vs eta scatter plot
  std::map<std::string, Histogram2> m_pullVsEta;
  /// Pull vs pT scatter plot
  std::map<std::string, Histogram2> m_pullVsPt;
  /// Pull vs eta-phi scatter plot
  std::map<std::string, Histogram3> m_pullVsEtaPhi;
  /// Pull vs eta-pT scatter plot
  std::map<std::string, Histogram3> m_pullVsEtaPt;
};

}  // namespace ActsExamples
