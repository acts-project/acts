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

static constexpr double nan = std::numeric_limits<double>::quiet_NaN();

ResPlotTool::ResPlotTool(const ResPlotTool::Config& cfg,
                         Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("ResPlotTool", lvl)) {
  const auto& etaAxis = m_cfg.varBinning.at("Eta");
  const auto& phiAxis = m_cfg.varBinning.at("Phi");
  const auto& ptAxis = m_cfg.varBinning.at("Pt");
  const auto& pullAxis = m_cfg.varBinning.at("Pull");

  ACTS_DEBUG("Initialize the histograms for residual and pull plots");

  std::vector<std::string> allParamNames = m_cfg.paramNames;
  allParamNames.push_back(m_cfg.qOverPtName);
  allParamNames.push_back(m_cfg.relQoverPtName);

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

    // residual vs eta-phi scatter plots
    m_resVsEtaPhi.emplace(parName,
                          Acts::Experimental::Histogram3(
                              std::format("res_{}_vs_eta_phi", parName),
                              std::format("Residual of {} vs eta-phi", parName),
                              std::array{etaAxis, phiAxis, residualAxis}));

    // residual vs eta-pT scatter plots
    m_resVsEtaPt.emplace(parName,
                         Acts::Experimental::Histogram3(
                             std::format("res_{}_vs_eta_pT", parName),
                             std::format("Residual of {} vs eta-pT", parName),
                             std::array{etaAxis, ptAxis, residualAxis}));

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

    // pull vs eta-phi scatter plots
    m_pullVsEtaPhi.emplace(parName,
                           Acts::Experimental::Histogram3(
                               std::format("pull_{}_vs_eta_phi", parName),
                               std::format("Pull of {} vs eta-phi", parName),
                               std::array{etaAxis, phiAxis, pullAxis}));

    // pull vs eta-pT scatter plots
    m_pullVsEtaPt.emplace(parName,
                          Acts::Experimental::Histogram3(
                              std::format("pull_{}_vs_eta_pT", parName),
                              std::format("Pull of {} vs eta-pT", parName),
                              std::array{etaAxis, ptAxis, pullAxis}));
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
  const ParametersVector& trackParameters = fittedParamters.parameters();
  const CovarianceMatrix& trackCovariance =
      fittedParamters.covariance().value_or(CovarianceMatrix::Zero());

  // get the perigee surface
  const Acts::Surface& pSurface = fittedParamters.referenceSurface();

  // get the truth parameter at the perigee surface
  ParametersVector truthParameters = ParametersVector::Zero();
  const Acts::Intersection3D intersection =
      pSurface
          .intersect(gctx, truthParticle.position(), truthParticle.direction())
          .closest();
  if (intersection.isValid()) {
    const Acts::Result<Acts::Vector2> lpResult = pSurface.globalToLocal(
        gctx, intersection.position(), truthParticle.direction());
    assert(lpResult.ok());

    truthParameters[eBoundLoc0] = lpResult.value()[eBoundLoc0];
    truthParameters[eBoundLoc1] = lpResult.value()[eBoundLoc1];
  } else {
    ACTS_ERROR("Cannot get the truth perigee parameter");
  }
  truthParameters[eBoundPhi] = phi(truthParticle.direction());
  truthParameters[eBoundTheta] = theta(truthParticle.direction());
  truthParameters[eBoundQOverP] = truthParticle.qOverP();
  truthParameters[eBoundTime] = truthParticle.time();

  // get the truth eta and pT
  const double truthEta = eta(truthParticle.direction());
  const double truthPhi = phi(truthParticle.direction());
  const double truthPt = truthParticle.transverseMomentum();

  // fill the histograms for residual and pull
  for (unsigned int paramId = 0; paramId < Acts::eBoundSize; paramId++) {
    const std::string& parName = m_cfg.paramNames.at(paramId);

    const double residual = trackParameters[paramId] - truthParameters[paramId];
    m_res.at(parName).fill({residual});
    m_resVsEta.at(parName).fill({truthEta, residual});
    m_resVsPt.at(parName).fill({truthPt, residual});
    m_resVsEtaPhi.at(parName).fill({truthEta, truthPhi, residual});
    m_resVsEtaPt.at(parName).fill({truthEta, truthPt, residual});

    const double var = trackCovariance(paramId, paramId);

    const double pull = var > 0 ? residual / std::sqrt(var) : nan;
    m_pull.at(parName).fill({pull});
    m_pullVsEta.at(parName).fill({truthEta, pull});
    m_pullVsPt.at(parName).fill({truthPt, pull});
    m_pullVsEtaPhi.at(parName).fill({truthEta, truthPhi, pull});
    m_pullVsEtaPt.at(parName).fill({truthEta, truthPt, pull});
  }

  // `reco(q/pT)` and `true(pT/q) * reco(q/pT)` residual and pull
  {
    const double truthQoverPt = truthParticle.charge() / truthPt;
    const double truthPtOverQ = truthPt / truthParticle.charge();
    const double recoQoverPt =
        trackParameters[eBoundQOverP] / std::sin(trackParameters[eBoundTheta]);
    const double residualQoverPt = recoQoverPt - truthQoverPt;
    m_res.at(m_cfg.qOverPtName).fill({residualQoverPt});
    m_resVsEta.at(m_cfg.qOverPtName).fill({truthEta, residualQoverPt});
    m_resVsPt.at(m_cfg.qOverPtName).fill({truthPt, residualQoverPt});

    const double residualRelQoverPt = truthPtOverQ * residualQoverPt;
    m_res.at(m_cfg.relQoverPtName).fill({residualRelQoverPt});
    m_resVsEta.at(m_cfg.relQoverPtName).fill({truthEta, residualRelQoverPt});
    m_resVsPt.at(m_cfg.relQoverPtName).fill({truthPt, residualRelQoverPt});

    const double covarianceQoverPt = [&]() {
      const Acts::Vector2 jacobian{
          -recoQoverPt / std::tan(trackParameters[eBoundTheta]),
          1 / std::sin(trackParameters[eBoundTheta])};
      const Acts::SquareMatrix2 covariance = trackCovariance(
          {eBoundTheta, eBoundQOverP}, {eBoundTheta, eBoundQOverP});
      return jacobian.transpose() * covariance * jacobian;
    }();
    const double covarianceRelQoverPt =
        Acts::square(truthPtOverQ) * covarianceQoverPt;

    const double pullQoverPt =
        covarianceQoverPt > 0 ? residualQoverPt / std::sqrt(covarianceQoverPt)
                              : nan;
    m_pull.at(m_cfg.qOverPtName).fill({pullQoverPt});
    m_pullVsEta.at(m_cfg.qOverPtName).fill({truthEta, pullQoverPt});
    m_pullVsPt.at(m_cfg.qOverPtName).fill({truthPt, pullQoverPt});

    const double pullRelQoverPt =
        covarianceRelQoverPt > 0
            ? residualRelQoverPt / std::sqrt(covarianceRelQoverPt)
            : nan;
    m_pull.at(m_cfg.relQoverPtName).fill({pullRelQoverPt});
    m_pullVsEta.at(m_cfg.relQoverPtName).fill({truthEta, pullRelQoverPt});
    m_pullVsPt.at(m_cfg.relQoverPtName).fill({truthPt, pullRelQoverPt});
  }
}

}  // namespace ActsExamples
