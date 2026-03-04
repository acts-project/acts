// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/ResPlotTool.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <format>

namespace ActsExamples {

ResPlotTool::ResPlotTool(const ResPlotTool::Config& cfg,
                         Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("ResPlotTool", lvl)) {
  const auto& etaAxis = m_cfg.varBinning.at("Eta");
  const auto& ptAxis = m_cfg.varBinning.at("Pt");
  const auto& pullAxis = m_cfg.varBinning.at("Pull");

  ACTS_DEBUG("Initialize the histograms for residual and pull plots");

  std::vector<std::string> allParamNames = m_cfg.paramNames;
  allParamNames.push_back(m_cfg.qOverPtName);
  allParamNames.push_back(m_cfg.ptQoverPtName);

  for (const std::string& parName : allParamNames) {
    const std::string parResidual = "Residual_" + parName;
    const auto& residualAxis = m_cfg.varBinning.at(parResidual);

    // residual distributions
    m_res.emplace(parName, Acts::Experimental::Histogram1(
                               std::format("res_{}", parName),
                               std::format("Residual of {}", parName),
                               std::array{residualAxis}));

    // residual vs eta scatter plots
    m_resVsEta.emplace(parName,
                       Acts::Experimental::Histogram2(
                           std::format("res_{}_vs_eta", parName),
                           std::format("Residual of {} vs eta", parName),
                           std::array{etaAxis, residualAxis}));

    // residual vs pT scatter plots
    m_resVsPt.emplace(parName, Acts::Experimental::Histogram2(
                                   std::format("res_{}_vs_pT", parName),
                                   std::format("Residual of {} vs pT", parName),
                                   std::array{ptAxis, residualAxis}));

    // pull distributions
    m_pull.emplace(
        parName, Acts::Experimental::Histogram1(
                     std::format("pull_{}", parName),
                     std::format("Pull of {}", parName), std::array{pullAxis}));

    // pull vs eta scatter plots
    m_pullVsEta.emplace(parName, Acts::Experimental::Histogram2(
                                     std::format("pull_{}_vs_eta", parName),
                                     std::format("Pull of {} vs eta", parName),
                                     std::array{etaAxis, pullAxis}));

    // pull vs pT scatter plots
    m_pullVsPt.emplace(parName, Acts::Experimental::Histogram2(
                                    std::format("pull_{}_vs_pT", parName),
                                    std::format("Pull of {} vs pT", parName),
                                    std::array{ptAxis, pullAxis}));
  }
}

void ResPlotTool::fill(const Acts::GeometryContext& gctx,
                       const SimParticleState& truthParticle,
                       const Acts::BoundTrackParameters& fittedParamters) {
  using ParametersVector = Acts::BoundTrackParameters::ParametersVector;
  using CovarianceMatrix = Acts::BoundTrackParameters::CovarianceMatrix;

  using Acts::VectorHelpers::eta;
  using Acts::VectorHelpers::perp;
  using Acts::VectorHelpers::phi;
  using Acts::VectorHelpers::theta;

  using enum Acts::BoundIndices;

  // get the fitted parameter (at perigee surface) and its error
  const ParametersVector& trackParameter = fittedParamters.parameters();
  const CovarianceMatrix& trackCovariance =
      fittedParamters.covariance().value_or(CovarianceMatrix::Zero());

  // get the perigee surface
  const Acts::Surface& pSurface = fittedParamters.referenceSurface();

  // get the truth parameter at the perigee surface
  ParametersVector truthParameter = ParametersVector::Zero();
  const Acts::Intersection3D intersection =
      pSurface
          .intersect(gctx, truthParticle.position(), truthParticle.direction())
          .closest();
  if (intersection.isValid()) {
    const Acts::Result<Acts::Vector2> lpResult = pSurface.globalToLocal(
        gctx, intersection.position(), truthParticle.direction());
    assert(lpResult.ok());

    truthParameter[eBoundLoc0] = lpResult.value()[eBoundLoc0];
    truthParameter[eBoundLoc1] = lpResult.value()[eBoundLoc1];
  } else {
    ACTS_ERROR("Cannot get the truth perigee parameter");
  }
  truthParameter[eBoundPhi] = phi(truthParticle.direction());
  truthParameter[eBoundTheta] = theta(truthParticle.direction());
  truthParameter[eBoundQOverP] = truthParticle.qOverP();
  truthParameter[eBoundTime] = truthParticle.time();

  // get the truth eta and pT
  const double truthEta = eta(truthParticle.direction());
  const double truthPt = truthParticle.transverseMomentum();

  // fill the histograms for residual and pull
  for (unsigned int parID = 0; parID < Acts::eBoundSize; parID++) {
    const std::string& parName = m_cfg.paramNames.at(parID);

    const double residual = trackParameter[parID] - truthParameter[parID];
    m_res.at(parName).fill({residual});
    m_resVsEta.at(parName).fill({truthEta, residual});
    m_resVsPt.at(parName).fill({truthPt, residual});

    if (!fittedParamters.covariance().has_value()) {
      ACTS_WARNING("Fitted track parameter :" << parName
                                              << " has no covariance");
      continue;
    }

    if (trackCovariance(parID, parID) <= 0.0) {
      ACTS_WARNING("Fitted track parameter :"
                   << parName << " has non-positive covariance = "
                   << trackCovariance(parID, parID));
      continue;
    }

    const double pull = residual / std::sqrt(trackCovariance(parID, parID));
    m_pull.at(parName).fill({pull});
    m_pullVsEta.at(parName).fill({truthEta, pull});
    m_pullVsPt.at(parName).fill({truthPt, pull});
  }

  // `q/pT` and `pT * q/pT` residual and pull
  do {
    const std::string& parName = m_cfg.ptQoverPtName;

    const double truthQoverPt = truthParticle.charge() / truthPt;
    const double recoQoverPt =
        trackParameter[eBoundQOverP] / std::sin(trackParameter[eBoundTheta]);
    const double residualQoverPt = recoQoverPt - truthQoverPt;
    m_res.at(m_cfg.qOverPtName).fill({residualQoverPt});
    m_resVsEta.at(m_cfg.qOverPtName).fill({truthEta, residualQoverPt});
    m_resVsPt.at(m_cfg.qOverPtName).fill({truthPt, residualQoverPt});

    const double residualPtQoverPt = truthQoverPt * residualQoverPt;
    m_res.at(m_cfg.ptQoverPtName).fill({residualPtQoverPt});
    m_resVsEta.at(m_cfg.ptQoverPtName).fill({truthEta, residualPtQoverPt});
    m_resVsPt.at(m_cfg.ptQoverPtName).fill({truthPt, residualPtQoverPt});

    const double covarianceQoverPt = [&]() {
      const Acts::Vector2 jacobian{
          -recoQoverPt / std::tan(trackParameter[eBoundTheta]),
          1 / std::sin(trackParameter[eBoundTheta])};
      const Acts::SquareMatrix2 covariance = trackCovariance(
          {eBoundTheta, eBoundQOverP}, {eBoundTheta, eBoundQOverP});
      return jacobian.transpose() * covariance * jacobian;
    }();
    const double covariancePtQoverPt =
        covarianceQoverPt * Acts::square(truthQoverPt);
    if (covariancePtQoverPt <= 0.0) {
      ACTS_WARNING("Fitted track parameter :"
                   << parName
                   << " has non-positive covariance = " << covariancePtQoverPt);
      continue;
    }

    const double pullQoverPt = residualQoverPt / std::sqrt(covarianceQoverPt);
    m_pull.at(m_cfg.qOverPtName).fill({pullQoverPt});
    m_pullVsEta.at(m_cfg.qOverPtName).fill({truthEta, pullQoverPt});
    m_pullVsPt.at(m_cfg.qOverPtName).fill({truthPt, pullQoverPt});

    const double pullPtQoverPt =
        residualPtQoverPt / std::sqrt(covariancePtQoverPt);
    m_pull.at(m_cfg.ptQoverPtName).fill({pullPtQoverPt});
    m_pullVsEta.at(m_cfg.ptQoverPtName).fill({truthEta, pullPtQoverPt});
    m_pullVsPt.at(m_cfg.ptQoverPtName).fill({truthPt, pullPtQoverPt});
  } while (false);
}

}  // namespace ActsExamples
