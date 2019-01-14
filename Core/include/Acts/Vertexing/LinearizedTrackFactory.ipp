// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Vertexing/LinearizedTrackFactory.hpp"

template <typename BField>
Acts::LinearizedTrack*
Acts::LinearizedTrackFactory<BField>::linearizeTrack(
    const Acts::BoundParameters* params,
    const Acts::Vector3D&        linPoint) const
{
  if (!params) return nullptr;

  const std::shared_ptr<Acts::PerigeeSurface> perigeeSurface
      = Acts::Surface::makeShared<Acts::PerigeeSurface>(linPoint);

  Acts::EigenStepper<BField> stepper(m_cfg.bField);

  // Set up propagator with void navigator
  Acts::Propagator<Acts::EigenStepper<BField>> propagator(stepper);

  // Set up propagator options
  Acts::PropagatorOptions<> options;
  // Variables to store track params and position at PCA to linPoint
  Acts::ActsVectorD<5>    paramsAtPCA;
  Acts::Vector3D          positionAtPCA;
  Acts::ActsSymMatrixD<5> parCovarianceAtPCA;

  // Do the propagation to linPoint
  const auto& result = propagator.propagate(*params, *perigeeSurface, options);
  if (result.status == Acts::Status::SUCCESS) {
    paramsAtPCA        = result.endParameters->parameters();
    positionAtPCA      = result.endParameters->position();
    parCovarianceAtPCA = *result.endParameters->covariance();
  } else {
    // *params are already perigeeParameters at linPoint (no propagation to
    // linPoint needed)
    paramsAtPCA        = params->parameters();
    positionAtPCA      = params->position();
    parCovarianceAtPCA = *params->covariance();
  }

  // phi_v and functions
  double phi_v     = paramsAtPCA(Acts::ParID_t::ePHI);
  double sin_phi_v = sin(phi_v);
  double cos_phi_v = cos(phi_v);

  // theta and functions
  double th     = paramsAtPCA(Acts::ParID_t::eTHETA);
  double sin_th = sin(th);
  double tan_th = tan(th);

  // q over p
  double q_ov_p = paramsAtPCA(Acts::ParID_t::eQOP);
  double sgn_h  = (q_ov_p < 0.) ? -1 : 1;

  Acts::Vector3D momentumAtPCA(phi_v, th, q_ov_p);

  // get B-field z-component at current position
  double B_z = m_cfg.bField.getField(linPoint)[Acts::eZ];

  double rho;
  // Curvature is infinite w/o b field
  if (B_z == 0. || std::abs(q_ov_p) < 1.e-15) {
    rho = 1.e+15;
  } else
    // signed(!) rho in _mm
    rho = sin_th * Acts::units::Nat2SI<Acts::units::MOMENTUM>
          (1 / q_ov_p * Acts::units::_GeV)
            / B_z * Acts::units::_mm; 

  double X  = positionAtPCA(0) - linPoint.x() + rho * sin_phi_v;
  double Y  = positionAtPCA(1) - linPoint.y() - rho * cos_phi_v;
  double S2 = (X * X + Y * Y);
  double S  = sqrt(S2);

  // F(V, p_i) at PCA (see Billoir paper)
  Acts::ActsVectorD<5> predParamsAtPCA;

  int sgnX = (X < 0.) ? -1 : 1;
  int sgnY = (Y < 0.) ? -1 : 1;

  double phiAtPCA;
  if (std::abs(X) > std::abs(Y))
    phiAtPCA = sgn_h * sgnX * std::acos(-sgn_h * Y / S);
  else {
    phiAtPCA                         = std::asin(sgn_h * X / S);
    if ((sgn_h * sgnY) > 0) phiAtPCA = sgn_h * sgnX * M_PI - phiAtPCA;
  }

  predParamsAtPCA[0] = rho - sgn_h * S;

  predParamsAtPCA[1] = positionAtPCA[Acts::eZ] - linPoint.z()
      + rho * (phi_v - phiAtPCA) / tan_th;
  predParamsAtPCA[2] = phiAtPCA;
  predParamsAtPCA[3] = th;
  predParamsAtPCA[4] = q_ov_p;

  // Fill position jacobian (D_k matrix)
  Acts::ActsMatrixD<5, 3> positionJacobian;
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

  // Fill momentum jacobian (E_k matrix)
  Acts::ActsMatrixD<5, 3> momentumJacobian;
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
  Acts::ActsVectorD<5> constTerm = predParamsAtPCA
      - positionJacobian * positionAtPCA - momentumJacobian * momentumAtPCA;

  return new Acts::LinearizedTrack(paramsAtPCA,
                                   parCovarianceAtPCA,
                                   linPoint,
                                   positionJacobian,
                                   momentumJacobian,
                                   positionAtPCA,
                                   momentumAtPCA,
                                   constTerm);
}
