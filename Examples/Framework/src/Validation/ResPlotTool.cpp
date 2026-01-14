// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/ResPlotTool.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <cmath>
#include <format>
#include <optional>
#include <ostream>

namespace ActsExamples {

ResPlotTool::ResPlotTool(const ResPlotTool::Config& cfg,
                         Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("ResPlotTool", lvl)) {}

void ResPlotTool::book(Cache& cache) const {
  const auto& etaAxis = m_cfg.varBinning.at("Eta");
  const auto& ptAxis = m_cfg.varBinning.at("Pt");
  const auto& pullAxis = m_cfg.varBinning.at("Pull");

  ACTS_DEBUG("Initialize the histograms for residual and pull plots");
  for (unsigned int parID = 0; parID < Acts::eBoundSize; parID++) {
    std::string parName = m_cfg.paramNames.at(parID);

    std::string parResidual = "Residual_" + parName;
    const auto& residualAxis = m_cfg.varBinning.at(parResidual);

    // residual distributions
    cache.res.emplace(parName, Acts::Experimental::Histogram1(
                                    std::format("res_{}", parName),
                                    std::format("Residual of {}", parName),
                                    std::array{residualAxis}));

    // residual vs eta scatter plots
    cache.res_vs_eta.emplace(
        parName,
        Acts::Experimental::Histogram2(std::format("res_{}_vs_eta", parName),
                                       std::format("Residual of {} vs eta", parName),
                                       std::array{etaAxis, residualAxis}));

    // residual vs pT scatter plots
    cache.res_vs_pT.emplace(
        parName,
        Acts::Experimental::Histogram2(std::format("res_{}_vs_pT", parName),
                                       std::format("Residual of {} vs pT", parName),
                                       std::array{ptAxis, residualAxis}));

    // pull distributions
    cache.pull.emplace(parName,
                       Acts::Experimental::Histogram1(
                           std::format("pull_{}", parName),
                           std::format("Pull of {}", parName),
                           std::array{pullAxis}));

    // pull vs eta scatter plots
    cache.pull_vs_eta.emplace(
        parName, Acts::Experimental::Histogram2(
                     std::format("pull_{}_vs_eta", parName),
                     std::format("Pull of {} vs eta", parName),
                     std::array{etaAxis, pullAxis}));

    // pull vs pT scatter plots
    cache.pull_vs_pT.emplace(
        parName, Acts::Experimental::Histogram2(
                     std::format("pull_{}_vs_pT", parName),
                     std::format("Pull of {} vs pT", parName),
                     std::array{ptAxis, pullAxis}));
  }
}

void ResPlotTool::fill(
    Cache& cache, const Acts::GeometryContext& gctx,
    const SimParticleState& truthParticle,
    const Acts::BoundTrackParameters& fittedParamters) const {
  using ParametersVector = Acts::BoundTrackParameters::ParametersVector;
  using Acts::VectorHelpers::eta;
  using Acts::VectorHelpers::perp;
  using Acts::VectorHelpers::phi;
  using Acts::VectorHelpers::theta;

  // get the fitted parameter (at perigee surface) and its error
  Acts::BoundVector trackParameter = fittedParamters.parameters();

  // get the perigee surface
  const Acts::Surface& pSurface = fittedParamters.referenceSurface();

  // get the truth position and momentum
  ParametersVector truthParameter = ParametersVector::Zero();

  // get the truth perigee parameter
  Acts::Intersection3D intersection =
      pSurface
          .intersect(gctx, truthParticle.position(), truthParticle.direction())
          .closest();
  if (intersection.isValid()) {
    auto lpResult = pSurface.globalToLocal(gctx, intersection.position(),
                                           truthParticle.direction());
    assert(lpResult.ok());

    truthParameter[Acts::BoundIndices::eBoundLoc0] =
        lpResult.value()[Acts::BoundIndices::eBoundLoc0];
    truthParameter[Acts::BoundIndices::eBoundLoc1] =
        lpResult.value()[Acts::BoundIndices::eBoundLoc1];
  } else {
    ACTS_ERROR("Cannot get the truth perigee parameter");
  }
  truthParameter[Acts::BoundIndices::eBoundPhi] =
      phi(truthParticle.direction());
  truthParameter[Acts::BoundIndices::eBoundTheta] =
      theta(truthParticle.direction());
  truthParameter[Acts::BoundIndices::eBoundQOverP] = truthParticle.qOverP();
  truthParameter[Acts::BoundIndices::eBoundTime] = truthParticle.time();

  // get the truth eta and pT
  const auto truthEta = eta(truthParticle.direction());
  const auto truthPt = truthParticle.transverseMomentum();

  // fill the histograms for residual and pull
  for (unsigned int parID = 0; parID < Acts::eBoundSize; parID++) {
    std::string parName = m_cfg.paramNames.at(parID);
    double residual = trackParameter[parID] - truthParameter[parID];
    cache.res.at(parName).fill({residual});
    cache.res_vs_eta.at(parName).fill({truthEta, residual});
    cache.res_vs_pT.at(parName).fill({truthPt, residual});

    if (fittedParamters.covariance().has_value()) {
      auto covariance = *fittedParamters.covariance();
      if (covariance(parID, parID) > 0) {
        double pull = residual / std::sqrt(covariance(parID, parID));
        cache.pull.at(parName).fill({pull});
        cache.pull_vs_eta.at(parName).fill({truthEta, pull});
        cache.pull_vs_pT.at(parName).fill({truthPt, pull});
      } else {
        ACTS_WARNING("Fitted track parameter :" << parName
                                                << " has negative covariance = "
                                                << covariance(parID, parID));
      }
    } else {
      ACTS_WARNING("Fitted track parameter :" << parName
                                              << " has no covariance");
    }
  }
}

}  // namespace ActsExamples
