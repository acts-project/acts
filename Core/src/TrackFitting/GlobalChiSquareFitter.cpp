// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GlobalChiSquareFitter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"

namespace Acts::Experimental {

// A projector used for scattering. By using Jacobian * phiThetaProjector one
// gets only the derivatives for the variables phi and theta.
const Acts::Matrix<eBoundSize, 2> phiThetaProjector = [] {
  Acts::Matrix<eBoundSize, 2> m = Acts::Matrix<eBoundSize, 2>::Zero();
  m(eBoundPhi, 0) = 1.0;
  m(eBoundTheta, 1) = 1.0;
  return m;
}();

// A projector used for scattering and energy loss. By using Jacobian *
// phiThetaQOverPProjector one gets the derivatives for the variables phi,
// theta, and q/p.
const Acts::Matrix<eBoundSize, 3> phiThetaQOverPProjector = [] {
  Acts::Matrix<eBoundSize, 3> m = Acts::Matrix<eBoundSize, 3>::Zero();
  m(eBoundPhi, 0) = 1.0;
  m(eBoundTheta, 1) = 1.0;
  m(eBoundQOverP, 2) = 1.0;
  return m;
}();

const Acts::Matrix<eBoundSize, 1> qOverPProjector = [] {
  Acts::Matrix<eBoundSize, 1> m = Acts::Matrix<eBoundSize, 1>::Zero();
  m(eBoundQOverP, 0) = 1.0;
  return m;
}();

void Gx2fMaterialProperties::updateParameters(
    const Eigen::VectorXd& deltaParamsExtended, const std::size_t stateIdx) {
  if (m_scatterer.isValid()) {
    m_scatPhi += deltaParamsExtended[stateIdx];
    m_scatTheta += deltaParamsExtended[stateIdx + 1ul];
  }
  if (m_eloss.isValid()) {
    const std::size_t eLossIdx = stateIdx + (nDim() - 1ul);
    m_lostEnergy += deltaParamsExtended[eLossIdx];
  }
}

void Gx2fMaterialProperties::updateTrackParameters(
    BoundTrackParameters& trackPars)  {
  if (m_scatterer.isValid()) {
    trackPars.parameters()[eBoundPhi] += deltaPhi();
    trackPars.parameters()[eBoundTheta] += deltaTheta();
  }
  if (m_eloss.isValid()) {
    const ParticleHypothesis& hypot = trackPars.particleHypothesis();
    trackPars.parameters()[eBoundQOverP] = hypot.qOverP(
        trackPars.absoluteMomentum() - m_lostEnergy, trackPars.charge());
    m_qOverPnew = trackPars.parameters()[eBoundQOverP];
  }
}

void Gx2fMaterialProperties::contributionToGx2fSums(
    Gx2fSystem& extendedSystem, const double theta, const std::size_t stateIdx,
    const Logger& logger) const {
  if (m_scatterer.isValid()) {
    const double invCovTheta = Acts::square(1. / m_scatterer.sigma());
    const double sinThetaLoc = std::sin(theta);
    const double invCovPhi = invCovTheta * Acts::square(sinThetaLoc);

    // Phi contribution
    extendedSystem.aMatrix()(stateIdx, stateIdx) += invCovPhi;
    extendedSystem.bVector()(stateIdx, 0) -= invCovPhi * deltaPhi();
    extendedSystem.chi2() += invCovPhi * Acts::square(deltaPhi());

    // Theta Contribution
    extendedSystem.aMatrix()(stateIdx + 1, stateIdx + 1) += invCovTheta;
    extendedSystem.bVector()(stateIdx + 1, 0) -= invCovTheta * deltaTheta();
    extendedSystem.chi2() += invCovTheta * Acts::square(deltaTheta());
    ACTS_VERBOSE("Scattering contributions in contributionToGx2fSums:\n"
                 << "     index: " << stateIdx << "-" << (stateIdx + 1) << " \n"
                 << "    invCov:        " << invCovPhi << "\n"
                 << "    sinThetaLoc:   " << sinThetaLoc << "\n"
                 << "    Phi:\n"
                 << "        scattering angle:     " << deltaPhi() << "\n"
                 << "        aMatrix contribution: " << invCovPhi << "\n"
                 << "        bVector contribution: " << invCovPhi * deltaPhi()
                 << "\n"
                 << "        chi2sum contribution: "
                 << invCovPhi * Acts::square(deltaPhi()) << "\n"
                 << "    Theta:\n"
                 << "        scattering angle:     " << deltaTheta() << "\n"
                 << "        aMatrix contribution: " << invCovTheta << "\n"
                 << "        bVector contribution: "
                 << invCovTheta * deltaTheta() << "\n"
                 << "        chi2sum contribution: "
                 << invCovTheta * Acts::square(deltaTheta()));
  }

  if (m_eloss.isValid()) {
    const std::size_t eLossIdx = stateIdx + (nDim() - 1ul);
    extendedSystem.aMatrix()(eLossIdx, eLossIdx) += Acts::square(1./m_eloss.lostSigma());
    extendedSystem.bVector()(eLossIdx, 0) -= Acts::square(1./m_eloss.lostSigma()) * (m_lostEnergy - m_eloss.lostEnergy());
    extendedSystem.chi2() += Acts::square(
        (m_lostEnergy - m_eloss.lostEnergy()) / m_eloss.lostSigma());
    ACTS_VERBOSE("Energy loss contributions in contributionToGx2fSums:\n"
                 << "       index: " << eLossIdx << " \n"
                 << "       measured loss: " << m_eloss.lostEnergy() << "+-"
                 << (m_eloss.lostSigma()) << ", \n"
                 << "       currently estimated loss: " << m_lostEnergy << "\n"
                 << "  --> chi2 contribution: "
                 << Acts::square((m_lostEnergy - m_eloss.lostEnergy()) /
                                 m_eloss.lostSigma()) << "\n"
                 << "  --> aMatrix contribution: "
                 << Acts::square(1./m_eloss.lostSigma()) << "\n"
                 << "  --> bVector contribution: "
                 << Acts::square(1./m_eloss.lostSigma()) *
                        (m_lostEnergy - m_eloss.lostEnergy()))
                        ;
  }
}

void updateGx2fParams(
    BoundTrackParameters& params, const Eigen::VectorXd& deltaParamsExtended,
    std::unordered_map<GeometryIdentifier, Gx2fMaterialProperties>& materialMap,
    const std::vector<GeometryIdentifier>& geoIdVector) {
  // update params
  params.parameters() +=
      deltaParamsExtended.topLeftCorner<eBoundSize, 1>().eval();

  std::size_t deltaPosition{eBoundSize};
  // update the scattering angles and energy loss.
  for (const GeometryIdentifier& geoId : geoIdVector) {
    const auto materialMapId = materialMap.find(geoId);
    assert(materialMapId != materialMap.end() &&
           "No material properties found for material surface.");

    materialMapId->second.updateParameters(deltaParamsExtended, deltaPosition);
    deltaPosition += materialMapId->second.nDim();
  }
}

void updateGx2fCovarianceParams(BoundMatrix& fullCovariancePredicted,
                                Gx2fSystem& extendedSystem) {
  // make invertible
  for (std::size_t i = 0; i < extendedSystem.nDims(); ++i) {
    if (extendedSystem.aMatrix()(i, i) == 0.) {
      extendedSystem.aMatrix()(i, i) = 1.;
    }
  }

  visit_measurement(extendedSystem.findRequiredNdf(), [&](auto N) {
    fullCovariancePredicted.topLeftCorner<N, N>() =
        extendedSystem.aMatrix().inverse().topLeftCorner<N, N>();
  });

  return;
}

void addMeasurementToGx2fSumsBackend(
    Gx2fSystem& extendedSystem,
    const std::vector<BoundMatrix>& jacobianFromStart,
    const std::vector<std::size_t>& materialIndices,
    const Eigen::MatrixXd& covarianceMeasurement, const BoundVector& predicted,
    const Eigen::VectorXd& measurement, const Eigen::MatrixXd& projector,
    const Logger& logger) {
  // First, we try to invert the covariance matrix. If the inversion fails, we
  // can already abort.
  const auto safeInvCovMeasurement = safeInverse(covarianceMeasurement);
  if (!safeInvCovMeasurement) {
    ACTS_WARNING("addMeasurementToGx2fSums: safeInvCovMeasurement failed.");
    ACTS_VERBOSE("    covarianceMeasurement:\n" << covarianceMeasurement);
    return;
  }
  ACTS_INFO("Jacobian size: " << jacobianFromStart.size());
  // Create an extended Jacobian. This one contains only eBoundSize rows,
  // because the rest is irrelevant. We fill it in the next steps.
  // TODO make dimsExtendedParams template with unrolling
  Eigen::MatrixXd extendedJacobian =
      Eigen::MatrixXd::Zero(eBoundSize, extendedSystem.nDims());

  // This part of the Jacobian comes from the material-less propagation
  extendedJacobian.topLeftCorner<eBoundSize, eBoundSize>() =
      jacobianFromStart[0];

  // If we have material, loop here over all Jacobians. We add extra columns for
  // their phi-theta projections. These parts account for the propagation of the
  // scattering angles. if energy loss is enabled, we also add the q/p
  // projection.
  for (std::size_t matSurface = 1; matSurface < jacobianFromStart.size();
       matSurface++) {
    const BoundMatrix& jac = jacobianFromStart[matSurface];

    // The position, where we need to insert the values in the extended Jacobian
    // accounting for material effects
    const std::size_t deltaPosition = materialIndices.at(matSurface - 1);
    ACTS_DEBUG("Update material jacobian " << deltaPosition << ", "
                                           << extendedSystem.nDims());

    //check the dimension of the material map: 2 means only scattering, 3 means scattering and energy loss
    const std::size_t matDim = materialIndices.at(matSurface) - deltaPosition;

    //check if energy loss is enabled, if yes, we need to add the q/p projection in the extended jacobian
    switch (matDim) {
      case 0:  // no material
        break;
      case 1:  // energy loss only
        extendedJacobian.block<eBoundSize, 1>(0, deltaPosition) =
            jac * qOverPProjector;
        break;
      case 2:  // scattering only
        extendedJacobian.block<eBoundSize, 2>(0, deltaPosition) =
            jac * phiThetaProjector;
        break;
      case 3:  // scattering and energy loss
        extendedJacobian.block<eBoundSize, 3>(0, deltaPosition) =
            jac * phiThetaQOverPProjector;
        break;
      default:
        ACTS_ERROR("Invalid material dimension: " << matDim<< " - this should be 1, 2, or 3.");
        throw std::domain_error("Invalid material dimension.");
    }
    
  }

  const Eigen::MatrixXd projJacobian = projector * extendedJacobian;

  const Eigen::VectorXd projPredicted = projector * predicted;

  const Eigen::VectorXd residual = measurement - projPredicted;

  // Finally contribute to chi2sum, aMatrix, and bVector
  extendedSystem.chi2() +=
      (residual.transpose() * (*safeInvCovMeasurement) * residual)(0, 0);

  extendedSystem.aMatrix() +=
      (projJacobian.transpose() * (*safeInvCovMeasurement) * projJacobian)
          .eval();

  extendedSystem.bVector() +=
      (residual.transpose() * (*safeInvCovMeasurement) * projJacobian)
          .eval()
          .transpose();

  ACTS_VERBOSE(
      "Contributions in addMeasurementToGx2fSums:\n"
      << "    predicted:   " << predicted.transpose() << "\n"
      << "    measurement: " << measurement.transpose() << "\n"
      << "    covarianceMeasurement:\n"
      << covarianceMeasurement << "\n"
      << "    projector:\n"
      << projector.eval() << "\n"
      << "    projJacobian:\n"
      << projJacobian.eval() << "\n"
      << "    projPredicted: " << (projPredicted.transpose()).eval() << "\n"
      << "    residual: " << (residual.transpose()).eval() << "\n"
      << "    extendedJacobian:\n"
      << extendedJacobian << "\n"
      << "    aMatrix contribution:\n"
      << (projJacobian.transpose() * (*safeInvCovMeasurement) * projJacobian)
             .eval()
      << "\n"
      << "    bVector contribution: "
      << (residual.transpose() * (*safeInvCovMeasurement) * projJacobian).eval()
      << "\n"
      << "    chi2sum contribution: "
      << (residual.transpose() * (*safeInvCovMeasurement) * residual)(0, 0)
      << "\n"
      << "    safeInvCovMeasurement:\n"
      << (*safeInvCovMeasurement));
} 

Eigen::VectorXd computeGx2fDeltaParams(const Gx2fSystem& extendedSystem) {
  return extendedSystem.aMatrix().colPivHouseholderQr().solve(
      extendedSystem.bVector());
}

}  // namespace Acts::Experimental
