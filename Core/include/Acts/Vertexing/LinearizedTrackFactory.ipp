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
  if (result.status == Status::SUCCESS) {
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

  // phi_v and functions
  double phi_v     = paramsAtPCA(ParID_t::ePHI);
  double sin_phi_v = sin(phi_v);
  double cos_phi_v = cos(phi_v);

  // theta and functions
  double th     = paramsAtPCA(ParID_t::eTHETA);
  double sin_th = sin(th);
  double tan_th = tan(th);

  // q over p
  double q_ov_p = paramsAtPCA(ParID_t::eQOP);
  double sgn_h  = (q_ov_p < 0.) ? -1 : 1;

  Vector3D momentumAtPCA(phi_v, th, q_ov_p);

  // get B-field z-component at current position
  double B_z = m_cfg.bField.getField(linPoint)[eZ];

  double rho;
  // Curvature is infinite w/o b field
  if (B_z == 0. || std::abs(q_ov_p) < 1.e-15) {
    rho = 1.e+15;
  } else {
    // signed(!) rho
    rho = sin_th * units::Nat2SI<units::MOMENTUM>(1 / q_ov_p) / B_z;
  }

  // Eq. 5.34 in Ref(1) (see .hpp)
  double X  = positionAtPCA(0) - linPoint.x() + rho * sin_phi_v;
  double Y  = positionAtPCA(1) - linPoint.y() - rho * cos_phi_v;
  double S2 = (X * X + Y * Y);
  double S  = sqrt(S2);

  /// F(V, p_i) at PCA in Billoir paper
  /// (see FullBilloirVertexFitter.hpp for paper reference,
  /// Page 140, Eq. (2) )
  ActsVectorD<5> predParamsAtPCA;

  int sgnX = (X < 0.) ? -1 : 1;
  int sgnY = (Y < 0.) ? -1 : 1;

  double phiAtPCA;
  if (std::abs(X) > std::abs(Y)) {
    phiAtPCA = sgn_h * sgnX * std::acos(-sgn_h * Y / S);
  } else {
    phiAtPCA = std::asin(sgn_h * X / S);
    if ((sgn_h * sgnY) > 0) {
      phiAtPCA = sgn_h * sgnX * M_PI - phiAtPCA;
    }
  }

  // Eq. 5.33 in Ref(1) (see .hpp)
  predParamsAtPCA[0] = rho - sgn_h * S;
  predParamsAtPCA[1]
      = positionAtPCA[eZ] - linPoint.z() + rho * (phi_v - phiAtPCA) / tan_th;
  predParamsAtPCA[2] = phiAtPCA;
  predParamsAtPCA[3] = th;
  predParamsAtPCA[4] = q_ov_p;

  // Fill position jacobian (D_k matrix), Eq. 5.36 in Ref(1)
  ActsMatrixD<5, 3> positionJacobian;
  positionJacobian.setZero();
  // First row
  positionJacobian(0, 0) = -sgn_h * X / S;
  positionJacobian(0, 1) = -sgn_h * Y / S;

  // Second row
  positionJacobian(1, 0) = rho * Y / (tan_th * S2);
  positionJacobian(1, 1) = -rho * X / (tan_th * S2);
  positionJacobian(1, 2) = 1.;

  // Third row
  positionJacobian(2, 0) = -Y / S2;
  positionJacobian(2, 1) = X / S2;

  // Fill momentum jacobian (E_k matrix), Eq. 5.37 in Ref(1)
  ActsMatrixD<5, 3> momentumJacobian;
  momentumJacobian.setZero();

  double R     = X * cos_phi_v + Y * sin_phi_v;
  double Q     = X * sin_phi_v - Y * cos_phi_v;
  double d_phi = phiAtPCA - phi_v;

  // First row
  momentumJacobian(0, 0) = -sgn_h * rho * R / S;

  double qOvS_red = 1 - sgn_h * Q / S;

  momentumJacobian(0, 1) = qOvS_red * rho / tan_th;
  momentumJacobian(0, 2) = -qOvS_red * rho / q_ov_p;

  // Second row
  momentumJacobian(1, 0) = (1 - rho * Q / S2) * rho / tan_th;
  momentumJacobian(1, 1) = (d_phi + rho * R / (S2 * tan_th * tan_th)) * rho;
  momentumJacobian(1, 2) = (d_phi - rho * R / S2) * rho / (q_ov_p * tan_th);

  // Third row
  momentumJacobian(2, 0) = rho * Q / S2;
  momentumJacobian(2, 1) = -rho * R / (S2 * tan_th);
  momentumJacobian(2, 2) = rho * R / (q_ov_p * S2);

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
