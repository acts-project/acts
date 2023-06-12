// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PerigeeSurface.hpp"

template <typename propagator_t, typename propagator_options_t>
Acts::Result<Acts::LinearizedTrack> Acts::
    HelicalTrackLinearizer<propagator_t, propagator_options_t>::linearizeTrack(
        const BoundTrackParameters& params, const Vector4& linPoint,
        const Acts::GeometryContext& gctx,
        const Acts::MagneticFieldContext& mctx, State& state) const {
  // Make Perigee surface at linPointPos, transverse plane of Perigee
  // corresponds the global x-y plane
  Vector3 linPointPos = VectorHelpers::position(linPoint);
  const std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(linPointPos);

  // Get intersection of the track with the Perigee if the particle would
  // move on a straight line.
  // This allows us to determine whether we need to propagate the track
  // forward or backward to arrive at the PCA.
  auto intersection = perigeeSurface->intersect(gctx, params.position(gctx),
                                                params.unitDirection(), false);

  // Create propagator options
  propagator_options_t pOptions(gctx, mctx);

  // Setting the propagation direction using the intersection length from
  // above
  // We handle zero path length as forward propagation, but we could actually
  // skip the whole propagation in this case
  pOptions.direction =
      Direction::fromScalarZeroAsPositive(intersection.intersection.pathLength);

  // Propagate to the PCA of linPointPos
  auto result = m_cfg.propagator->propagate(params, *perigeeSurface, pOptions);
  if (not result.ok()) {
    return result.error();
  }

  // Extracting the track parameters at said PCA - this corresponds to the
  // Perigee representation of the track wrt linPointPos
  const auto& endParams = *result->endParameters;
  BoundVector paramsAtPCA = endParams.parameters();

  // Extracting the 4D position of the PCA in global coordinates
  Vector4 pca = Vector4::Zero();
  {
    auto pos = endParams.position(gctx);
    pca[ePos0] = pos[ePos0];
    pca[ePos1] = pos[ePos1];
    pca[ePos2] = pos[ePos2];
    pca[eTime] = endParams.time();
  }
  BoundSymMatrix parCovarianceAtPCA = endParams.covariance().value();

  // Extracting Perigee parameters and compute functions of them for later
  // usage
  ActsScalar d0 = paramsAtPCA(BoundIndices::eBoundLoc0);

  ActsScalar phi = paramsAtPCA(BoundIndices::eBoundPhi);
  ActsScalar sinPhi = std::sin(phi);
  ActsScalar cosPhi = std::cos(phi);

  ActsScalar theta = paramsAtPCA(BoundIndices::eBoundTheta);
  ActsScalar sinTheta = std::sin(theta);
  ActsScalar tanTheta = std::tan(theta);

  // q over p
  ActsScalar qOvP = paramsAtPCA(BoundIndices::eBoundQOverP);

  // Get mass hypothesis from propagator options
  ActsScalar m0 = pOptions.mass;
  // Assume unit charge
  ActsScalar p = std::abs(1 / qOvP);
  // Speed in units of c
  ActsScalar beta = p / std::hypot(p, m0);
  // Transverse speed (i.e., speed in the x-y plane)
  ActsScalar betaT = beta * sinTheta;

  // Momentu at the PCA
  Vector3 momentumAtPCA(phi, theta, qOvP);

  // Define Jacobians, which will be filled later
  ActsMatrix<eBoundSize, 4> positionJacobian;
  positionJacobian.setZero();
  ActsMatrix<eBoundSize, 3> momentumJacobian;
  momentumJacobian.setZero();
  // get the z-component of the B-field at the PCA
  auto field =
      m_cfg.bField->getField(VectorHelpers::position(pca), state.fieldCache);
  if (!field.ok()) {
    return field.error();
  }
  ActsScalar Bz = (*field)[eZ];

  // If there is no magnetic field the particle has a straight trajectory
  // If there is a constant magnetic field it has a helical trajectory
  if (Bz == 0. || std::abs(qOvP) < m_cfg.minQoP) {
    // Fill position Jacobian, i.e., matrix A from Eq. 5.39 in Ref(1)
    // First row
    positionJacobian(0, 0) = -sinPhi;
    positionJacobian(0, 1) = cosPhi;

    // Second row
    positionJacobian(1, 0) = -cosPhi / tanTheta;
    positionJacobian(1, 1) = -sinPhi / tanTheta;
    positionJacobian(1, 2) = 1.;

    // Sixth row
    positionJacobian(5, 0) = -cosPhi / betaT;
    positionJacobian(5, 1) = -sinPhi / betaT;
    positionJacobian(5, 3) = 1.;

    // Fill momentum Jacobian, i.e., matrix B from Eq. 5.40 in Ref(1)
    // Second row
    momentumJacobian(1, 0) = -d0 / tanTheta;

    // Third row
    momentumJacobian(2, 0) = 1.;

    // Fourth row
    momentumJacobian(3, 1) = 1.;

    // Fifth row
    momentumJacobian(4, 2) = 1.;

    // Sixth row
    momentumJacobian(5, 0) = -d0 / betaT;
  } else {
    // Helix radius
    ActsScalar rho = sinTheta * (1. / qOvP) / Bz;
    // Sign of helix radius
    ActsScalar h = (rho < 0.) ? -1 : 1;

    // Quantities from Eq. 5.34 in Ref(1) (see .hpp)
    ActsScalar X = pca(0) - linPointPos.x() + rho * sinPhi;
    ActsScalar Y = pca(1) - linPointPos.y() - rho * cosPhi;
    ActsScalar S2 = (X * X + Y * Y);
    // S is the 2D distance from the helix center to linPointPos
    // in the x-y plane
    ActsScalar S = std::sqrt(S2);

    // Fill position Jacobian, i.e., matrix A from Eq. 5.36 in Ref(1)
    // First row
    positionJacobian(0, 0) = -h * X / S;
    positionJacobian(0, 1) = -h * Y / S;

    ActsScalar XoverS2 = X / S2;
    ActsScalar YoverS2 = Y / S2;
    ActsScalar rhoCotTheta = rho / tanTheta;
    ActsScalar rhoOverBetaT = rho / betaT;

    // Second row
    positionJacobian(1, 0) = rhoCotTheta * YoverS2;
    positionJacobian(1, 1) = -rhoCotTheta * XoverS2;
    positionJacobian(1, 2) = 1.;

    // Third row
    positionJacobian(2, 0) = -YoverS2;
    positionJacobian(2, 1) = XoverS2;

    // Sixth row
    positionJacobian(5, 0) = rhoOverBetaT * YoverS2;
    positionJacobian(5, 1) = -rhoOverBetaT * XoverS2;
    positionJacobian(5, 3) = 1.;

    // Fill momentum Jacobian, i.e., B matrix from Eq. 5.36 in Ref(1).
    // Since we propagated to the PCA (point P in Ref(1)), the points
    // P and V coincide and we can choose deltaPhi = 0.
    // One can show that if deltaPhi = 0 --> R = 0 and Q = h * S.
    // As a consequence, many terms of the B matrix from Eq. 5.36 vanish.

    // Absolute value of rho over S
    ActsScalar absRhoOverS = h * rho / S;

    // Second row
    momentumJacobian(1, 0) = rhoCotTheta * (1. - absRhoOverS);

    // Third row
    momentumJacobian(2, 0) = absRhoOverS;

    // Fourth and fifth row
    momentumJacobian(3, 1) = 1.;
    momentumJacobian(4, 2) = 1.;

    // Sixth row
    momentumJacobian(5, 0) = rhoOverBetaT * (1. - absRhoOverS);
  }

  // const term in Taylor expansion from Eq. 5.38 in Ref(1)
  BoundVector constTerm =
      paramsAtPCA - positionJacobian * pca - momentumJacobian * momentumAtPCA;

  // The parameter weight
  BoundSymMatrix weightAtPCA = parCovarianceAtPCA.inverse();

  return LinearizedTrack(paramsAtPCA, parCovarianceAtPCA, weightAtPCA, linPoint,
                         positionJacobian, momentumJacobian, pca, momentumAtPCA,
                         constTerm);
}
