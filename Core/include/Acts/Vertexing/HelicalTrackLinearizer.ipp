// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PerigeeSurface.hpp"

template <typename propagator_t, typename propagator_options_t>
Acts::Result<Acts::LinearizedTrack> Acts::
    HelicalTrackLinearizer<propagator_t, propagator_options_t>::linearizeTrack(
        const BoundParameters& params, const SpacePointVector& linPoint) const {
  Vector3D linPointPos = VectorHelpers::position(linPoint);

  const std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(linPointPos);

  const BoundParameters* endParams = nullptr;
  // Do the propagation to linPointPos
  auto result =
      m_cfg.propagator->propagate(params, *perigeeSurface, m_cfg.pOptions);
  if (result.ok()) {
    endParams = (*result).endParameters.get();

  } else {
    return result.error();
  }

  BoundVector paramsAtPCA = endParams->parameters();
  SpacePointVector positionAtPCA = SpacePointVector::Zero();
  VectorHelpers::position(positionAtPCA) = endParams->position();
  BoundSymMatrix parCovarianceAtPCA = *(endParams->covariance());

  if (endParams->covariance()->determinant() == 0) {
    // Use the original parameters
    paramsAtPCA = params.parameters();
    VectorHelpers::position(positionAtPCA) = params.position();
    parCovarianceAtPCA = *(params.covariance());
  }

  // phiV and functions
  double phiV = paramsAtPCA(ParID_t::ePHI);
  double sinPhiV = std::sin(phiV);
  double cosPhiV = std::cos(phiV);

  // theta and functions
  double th = paramsAtPCA(ParID_t::eTHETA);
  const double sinTh = std::sin(th);
  const double tanTh = std::tan(th);

  // q over p
  double qOvP = paramsAtPCA(ParID_t::eQOP);
  double sgnH = (qOvP < 0.) ? -1 : 1;

  Vector3D momentumAtPCA(phiV, th, qOvP);

  // get B-field z-component at current position
  double Bz = m_cfg.bField.getField(VectorHelpers::position(positionAtPCA))[eZ];

  double rho;
  double rho2;
  // Curvature is infinite w/o b field
  if (Bz == 0. || std::abs(qOvP) < m_cfg.minQoP) {
    rho = m_cfg.maxRho;
  } else {
    // signed(!) rho
    rho2 = sinTh / (qOvP * Bz);
    rho = sinTh * (1. / qOvP) / Bz;
  }

  // Eq. 5.34 in Ref(1) (see .hpp)
  double X = positionAtPCA(0) - linPointPos.x() + rho * sinPhiV;
  double Y = positionAtPCA(1) - linPointPos.y() - rho * cosPhiV;
  const double S2 = (X * X + Y * Y);
  const double S = std::sqrt(S2);

  /// F(V, p_i) at PCA in Billoir paper
  /// (see FullBilloirVertexFitter.hpp for paper reference,
  /// Page 140, Eq. (2) )
  BoundVector predParamsAtPCA;

  int sgnX = (X < 0.) ? -1 : 1;
  int sgnY = (Y < 0.) ? -1 : 1;

  double phiAtPCA;
  if (std::abs(X) > std::abs(Y)) {
    phiAtPCA = sgnH * sgnX * std::acos(-sgnH * Y / S);
  } else {
    phiAtPCA = std::asin(sgnH * X / S);
    if ((sgnH * sgnY) > 0) {
      phiAtPCA = sgnH * sgnX * M_PI - phiAtPCA;
    }
  }

  // Eq. 5.33 in Ref(1) (see .hpp)
  predParamsAtPCA[0] = rho - sgnH * S;
  predParamsAtPCA[1] =
      positionAtPCA[eZ] - linPointPos.z() + rho * (phiV - phiAtPCA) / tanTh;
  predParamsAtPCA[2] = phiAtPCA;
  predParamsAtPCA[3] = th;
  predParamsAtPCA[4] = qOvP;
  predParamsAtPCA[5] = 0.;

  // Fill position jacobian (D_k matrix), Eq. 5.36 in Ref(1)
  SpacePointToBoundMatrix positionJacobian;
  positionJacobian.setZero();
  // First row
  positionJacobian(0, 0) = -sgnH * X / S;
  positionJacobian(0, 1) = -sgnH * Y / S;

  const double S2tanTh = S2 * tanTh;

  // Second row
  positionJacobian(1, 0) = rho * Y / S2tanTh;
  positionJacobian(1, 1) = -rho * X / S2tanTh;
  positionJacobian(1, 2) = 1.;

  // Third row
  positionJacobian(2, 0) = -Y / S2;
  positionJacobian(2, 1) = X / S2;

  // TODO: include timing in track linearization
  // Last row
  positionJacobian(5, 3) = 1;

  // Fill momentum jacobian (E_k matrix), Eq. 5.37 in Ref(1)
  ActsMatrixD<BoundParsDim, 3> momentumJacobian;
  momentumJacobian.setZero();

  double R = X * cosPhiV + Y * sinPhiV;
  double Q = X * sinPhiV - Y * cosPhiV;
  double dPhi = phiAtPCA - phiV;

  // First row
  momentumJacobian(0, 0) = -sgnH * rho * R / S;

  double qOvSred = 1 - sgnH * Q / S;

  momentumJacobian(0, 1) = qOvSred * rho / tanTh;
  momentumJacobian(0, 2) = -qOvSred * rho / qOvP;

  const double rhoOverS2 = rho / S2;

  // Second row
  momentumJacobian(1, 0) = (1 - rhoOverS2 * Q) * rho / tanTh;
  momentumJacobian(1, 1) = (dPhi + rho * R / (S2tanTh * tanTh)) * rho;
  momentumJacobian(1, 2) = (dPhi - rhoOverS2 * R) * rho / (qOvP * tanTh);

  // Third row
  momentumJacobian(2, 0) = rhoOverS2 * Q;
  momentumJacobian(2, 1) = -rho * R / S2tanTh;
  momentumJacobian(2, 2) = rhoOverS2 * R / qOvP;

  // Last two rows:
  momentumJacobian(3, 1) = 1.;
  momentumJacobian(4, 2) = 1.;

  // const term F(V_0, p_0) in Talyor expansion
  BoundVector constTerm = predParamsAtPCA - positionJacobian * positionAtPCA -
                          momentumJacobian * momentumAtPCA;

  // The parameter weight
  ActsSymMatrixD<5> parWeight =
      (parCovarianceAtPCA.block<5, 5>(0, 0)).inverse();

  BoundSymMatrix weightAtPCA{BoundSymMatrix::Identity()};
  weightAtPCA.block<5, 5>(0, 0) = parWeight;

  return LinearizedTrack(paramsAtPCA, parCovarianceAtPCA, weightAtPCA, linPoint,
                         positionJacobian, momentumJacobian, positionAtPCA,
                         momentumAtPCA, constTerm);
}
