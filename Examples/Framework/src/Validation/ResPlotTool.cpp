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
#include <cstdint>
#include <format>
#include <stdexcept>

namespace ActsExamples {

static constexpr double nan = std::numeric_limits<double>::quiet_NaN();

ResPlotTool::ResPlotTool(const ResPlotTool::Config& cfg,
                         Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("ResPlotTool", lvl)) {
  // `varBinning.at` would only report the key type, not the missing key
  const auto binning = [this](const std::string& key) -> const AxisVariant& {
    const auto it = m_cfg.varBinning.find(key);
    if (it == m_cfg.varBinning.end()) {
      throw std::invalid_argument("ResPlotTool: missing binning for '" + key +
                                  "'");
    }
    return it->second;
  };

  const auto& etaAxis = binning("Eta");
  const auto& phiAxis = binning("Phi");
  const auto& ptAxis = binning("Pt");
  const auto& pullAxis = binning("Pull");

  if (m_cfg.paramNames.size() != Acts::eBoundSize) {
    throw std::invalid_argument(
        "ResPlotTool: expected one name per bound parameter");
  }

  ACTS_DEBUG("Initialize the histograms for residual and pull plots");

  std::vector<std::string> allParamNames = m_cfg.paramNames;
  allParamNames.push_back(m_cfg.qOverPtName);
  allParamNames.push_back(m_cfg.relQoverPtName);

  for (const std::string& parName : allParamNames) {
    const auto& residualAxis = binning("Residual_" + parName);

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
  using Acts::VectorHelpers::eta;
  using Acts::VectorHelpers::perp;
  using Acts::VectorHelpers::phi;
  using Acts::VectorHelpers::theta;

  using enum Acts::BoundIndices;

  // get the perigee surface
  const Acts::Surface& pSurface = fittedParamters.referenceSurface();

  // get the truth parameter at the perigee surface
  Acts::BoundVector truthParameters = Acts::BoundVector::Zero();
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

  // bin on the particle, not the bound parameters, which would round-trip the
  // direction through phi/theta
  fill(truthParameters,
       Binning{eta(truthParticle.direction()), phi(truthParticle.direction()),
               truthParticle.transverseMomentum()},
       truthParticle.charge(), truthParticle.absoluteCharge(), fittedParamters);
}

void ResPlotTool::fill(const Acts::BoundTrackParameters& truthParameters,
                       const Acts::BoundTrackParameters& fittedParameters) {
  using Acts::VectorHelpers::eta;
  using Acts::VectorHelpers::phi;

  if (truthParameters.referenceSurface() !=
      fittedParameters.referenceSurface()) {
    throw std::invalid_argument(
        "ResPlotTool: truth and fitted parameters are expressed on different "
        "reference surfaces");
  }

  const double truthCharge = truthParameters.charge();
  fill(truthParameters.parameters(),
       Binning{eta(truthParameters.direction()),
               phi(truthParameters.direction()),
               truthParameters.transverseMomentum()},
       truthCharge, std::abs(truthCharge), fittedParameters);
}

void ResPlotTool::fill(const Acts::BoundVector& truthVector,
                       const Binning& binning, double truthCharge,
                       double truthAbsCharge,
                       const Acts::BoundTrackParameters& fittedParameters) {
  using enum Acts::BoundIndices;

  const double truthEta = binning.eta;
  const double truthPhi = binning.phi;
  const double truthPt = binning.pt;

  // get the fitted parameter and its error
  const Acts::BoundVector& trackParameters = fittedParameters.parameters();
  const Acts::BoundMatrix& trackCovariance =
      fittedParameters.covariance().value_or(Acts::BoundMatrix::Zero());

  // fill the histograms for residual and pull
  for (unsigned int paramId = 0; paramId < Acts::eBoundSize; paramId++) {
    const std::string& parName = m_cfg.paramNames.at(paramId);

    const double residual = trackParameters[paramId] - truthVector[paramId];
    fillResidual(parName, residual, truthEta, truthPhi, truthPt);

    const double var = trackCovariance(paramId, paramId);

    const double pull = var > 0 ? residual / std::sqrt(var) : nan;
    fillPull(parName, pull, truthEta, truthPhi, truthPt);
  }

  // `reco(q/pT)` and `true(pT/q) * reco(q/pT)` residual and pull
  {
    const double truthQoverPt = truthCharge / truthPt;
    const double truthPtOverAbsQ = truthPt / truthAbsCharge;
    const double recoQoverPt =
        trackParameters[eBoundQOverP] / std::sin(trackParameters[eBoundTheta]);
    const double residualQoverPt = recoQoverPt - truthQoverPt;
    fillResidual(m_cfg.qOverPtName, residualQoverPt, truthEta, truthPhi,
                 truthPt);

    const double residualRelQoverPt = truthPtOverAbsQ * residualQoverPt;
    fillResidual(m_cfg.relQoverPtName, residualRelQoverPt, truthEta, truthPhi,
                 truthPt);

    const double covarianceQoverPt = [&]() -> double {
      const Acts::Vector2 jacobian{
          -recoQoverPt / std::tan(trackParameters[eBoundTheta]),
          1 / std::sin(trackParameters[eBoundTheta])};
      const Acts::SquareMatrix2 covariance = trackCovariance(
          {eBoundTheta, eBoundQOverP}, {eBoundTheta, eBoundQOverP});
      return jacobian.transpose() * covariance * jacobian;
    }();
    const double covarianceRelQoverPt =
        Acts::square(truthPtOverAbsQ) * covarianceQoverPt;

    const double pullQoverPt =
        covarianceQoverPt > 0 ? residualQoverPt / std::sqrt(covarianceQoverPt)
                              : nan;
    fillPull(m_cfg.qOverPtName, pullQoverPt, truthEta, truthPhi, truthPt);

    const double pullRelQoverPt =
        covarianceRelQoverPt > 0
            ? residualRelQoverPt / std::sqrt(covarianceRelQoverPt)
            : nan;
    fillPull(m_cfg.relQoverPtName, pullRelQoverPt, truthEta, truthPhi, truthPt);
  }
}

void ResPlotTool::fill(const Binning& binning,
                       const Acts::VariableBoundSubspaceHelper& subspace,
                       const Acts::BoundVector& residuals,
                       const Acts::BoundMatrix& residualCovariance) {
  for (const std::uint8_t index : subspace) {
    const std::string& parName = m_cfg.paramNames.at(index);

    const double residual = residuals[index];
    fillResidual(parName, residual, binning.eta, binning.phi, binning.pt);

    // `V - HPH^T` is not positive definite, so there is not always a pull
    const double var = residualCovariance(index, index);
    const double pull = var > 0 ? residual / std::sqrt(var) : nan;
    fillPull(parName, pull, binning.eta, binning.phi, binning.pt);
  }
}

void ResPlotTool::fillResidual(const std::string& paramName, double residual,
                               double truthEta, double truthPhi,
                               double truthPt) {
  m_res.at(paramName).fill({residual});
  m_resVsEta.at(paramName).fill({truthEta, residual});
  m_resVsPt.at(paramName).fill({truthPt, residual});
  m_resVsEtaPhi.at(paramName).fill({truthEta, truthPhi, residual});
  m_resVsEtaPt.at(paramName).fill({truthEta, truthPt, residual});
}

void ResPlotTool::fillPull(const std::string& paramName, double pull,
                           double truthEta, double truthPhi, double truthPt) {
  m_pull.at(paramName).fill({pull});
  m_pullVsEta.at(paramName).fill({truthEta, pull});
  m_pullVsPt.at(paramName).fill({truthPt, pull});
  m_pullVsEtaPhi.at(paramName).fill({truthEta, truthPhi, pull});
  m_pullVsEtaPt.at(paramName).fill({truthEta, truthPt, pull});
}

}  // namespace ActsExamples
