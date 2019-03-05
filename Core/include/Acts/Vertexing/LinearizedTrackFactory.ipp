// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PerigeeSurface.hpp"

template <typename BField,
          typename Propagator_t,
          typename action_list_t,
          typename aborter_list_t>
Acts::LinearizedTrack
Acts::LinearizedTrackFactory<BField,
                             Propagator_t,
                             action_list_t,
                             aborter_list_t>::
    linearizeTrack(const BoundParameters* params,
                   const Vector3D&        linPoint,
                   const Propagator_t&    propagator) const
{
  if (params == nullptr) {
    return LinearizedTrack();
  }

  const std::shared_ptr<PerigeeSurface> perigeeSurface
      = Surface::makeShared<PerigeeSurface>(linPoint);

  // Variables to store track params and position at PCA to linPoint
  ActsVectorD<5>    paramsAtPCA;
  Vector3D          positionAtPCA;
  ActsSymMatrixD<5> parCovarianceAtPCA;

  // Do the propagation to linPoint
  const auto& result
      = propagator.propagate(*params, *perigeeSurface, m_cfg.propagatorOptions);
  if (result.status == PropagatorStatus::SUCCESS) {
    paramsAtPCA        = result.endParameters->parameters();
    positionAtPCA      = result.endParameters->position();
    parCovarianceAtPCA = *result.endParameters->covariance();

    // *params are already perigeeParameters at linPoint (no propagation to
    // linPoint needed)
  } else {
    paramsAtPCA        = params->parameters();
    positionAtPCA      = params->position();
    parCovarianceAtPCA = *params->covariance();
  }

  // phiV and functions
  double phiV    = paramsAtPCA(ParID_t::ePHI);
  double sinPhiV = std::sin(phiV);
  double cosPhiV = std::cos(phiV);

  // theta and functions
  double th    = paramsAtPCA(ParID_t::eTHETA);
  double sinTh = std::sin(th);
  double tanTh = std::tan(th);

  // q over p
  double qOvP = paramsAtPCA(ParID_t::eQOP);
  double sgnH = (qOvP < 0.) ? -1 : 1;

  Vector3D momentumAtPCA(phiV, th, qOvP);

  // get B-field z-component at current position
  double Bz = m_cfg.bField.getField(linPoint)[eZ];

  double rho;
  // Curvature is infinite w/o b field
  if (Bz == 0. || std::abs(qOvP) < 1.e-15) {
    rho = std::numeric_limits<double>::max();
  } else {
    // signed(!) rho
    rho = sinTh * units::Nat2SI<units::MOMENTUM>(1 / qOvP) / Bz;
  }

  // Eq. 5.34 in Ref(1) (see .hpp)
  double X  = positionAtPCA(0) - linPoint.x() + rho * sinPhiV;
  double Y  = positionAtPCA(1) - linPoint.y() - rho * cosPhiV;
  double S2 = (X * X + Y * Y);
  double S  = std::sqrt(S2);

  /// F(V, p_i) at PCA in Billoir paper
  /// (see FullBilloirVertexFitter.hpp for paper reference,
  /// Page 140, Eq. (2) )
  ActsVectorD<5> predParamsAtPCA;

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
  predParamsAtPCA[1]
      = positionAtPCA[eZ] - linPoint.z() + rho * (phiV - phiAtPCA) / tanTh;
  predParamsAtPCA[2] = phiAtPCA;
  predParamsAtPCA[3] = th;
  predParamsAtPCA[4] = qOvP;

  // Fill position jacobian (D_k matrix), Eq. 5.36 in Ref(1)
  ActsMatrixD<5, 3> positionJacobian;
  positionJacobian.setZero();
  // First row
  positionJacobian(0, 0) = -sgnH * X / S;
  positionJacobian(0, 1) = -sgnH * Y / S;

  // Second row
  positionJacobian(1, 0) = rho * Y / (tanTh * S2);
  positionJacobian(1, 1) = -rho * X / (tanTh * S2);
  positionJacobian(1, 2) = 1.;

  // Third row
  positionJacobian(2, 0) = -Y / S2;
  positionJacobian(2, 1) = X / S2;

  // Fill momentum jacobian (E_k matrix), Eq. 5.37 in Ref(1)
  ActsMatrixD<5, 3> momentumJacobian;
  momentumJacobian.setZero();

  double R    = X * cosPhiV + Y * sinPhiV;
  double Q    = X * sinPhiV - Y * cosPhiV;
  double dPhi = phiAtPCA - phiV;

  // First row
  momentumJacobian(0, 0) = -sgnH * rho * R / S;

  double qOvSred = 1 - sgnH * Q / S;

  momentumJacobian(0, 1) = qOvSred * rho / tanTh;
  momentumJacobian(0, 2) = -qOvSred * rho / qOvP;

  // Second row
  momentumJacobian(1, 0) = (1 - rho * Q / S2) * rho / tanTh;
  momentumJacobian(1, 1) = (dPhi + rho * R / (S2 * tanTh * tanTh)) * rho;
  momentumJacobian(1, 2) = (dPhi - rho * R / S2) * rho / (qOvP * tanTh);

  // Third row
  momentumJacobian(2, 0) = rho * Q / S2;
  momentumJacobian(2, 1) = -rho * R / (S2 * tanTh);
  momentumJacobian(2, 2) = rho * R / (qOvP * S2);

  // Last two rows:
  momentumJacobian(3, 1) = 1.;
  momentumJacobian(4, 2) = 1.;

  // const term F(V_0, p_0) in Talyor expansion
  ActsVectorD<5> constTerm = predParamsAtPCA - positionJacobian * positionAtPCA
      - momentumJacobian * momentumAtPCA;

  return LinearizedTrack(paramsAtPCA,
                         parCovarianceAtPCA,
                         linPoint,
                         positionJacobian,
                         momentumJacobian,
                         positionAtPCA,
                         momentumAtPCA,
                         constTerm);
}
