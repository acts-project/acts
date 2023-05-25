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
  double phi = paramsAtPCA(BoundIndices::eBoundPhi);
  double sinPhi = std::sin(phi);
  double cosPhi = std::cos(phi);

  double theta = paramsAtPCA(BoundIndices::eBoundTheta);
  const double sinTheta = std::sin(theta);
  const double tanTheta = std::tan(theta);

  // q over p
  double qOvP = paramsAtPCA(BoundIndices::eBoundQOverP);

  // Charge of the particle, determines the sign of the helix radius rho
  double sgnH = (qOvP < 0.) ? -1 : 1;

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
  double Bz = (*field)[eZ];

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

    // TODO: include timing in track linearization - will be added
    // in next PR
    // Sixth row
    positionJacobian(5, 3) = 1.;

    // Quantities from Eq. 5.41 and 5.42 in Ref (1)
    double R = (pca(0) - linPointPos.x()) * cosPhi +
               (pca(1) - linPointPos.y()) * sinPhi;
    double Q = (pca(0) - linPointPos.x()) * sinPhi -
               (pca(1) - linPointPos.y()) * cosPhi;

    // Fill momentum Jacobian, i.e., matrix B from Eq. 5.40 in Ref(1)
    // First row
    momentumJacobian(0, 0) = -R;

    // Second row
    momentumJacobian(1, 0) = Q / tanTheta;
    momentumJacobian(1, 1) = R / (sinTheta * sinTheta);

    // Third row
    momentumJacobian(2, 0) = 1.;

    // Fourth row
    momentumJacobian(3, 1) = 1.;

    // Fifth row
    momentumJacobian(4, 2) = 1.;

    // TODO: Derivatives of time --> Next PR
  } else {
    // Helix radius
    double rho{sinTheta * (1. / qOvP) / Bz};

    // Quantities from Eq. 5.34 in Ref(1) (see .hpp)
    double X = pca(0) - linPointPos.x() + rho * sinPhi;
    double Y = pca(1) - linPointPos.y() - rho * cosPhi;
    const double S2 = (X * X + Y * Y);
    // S is the 2D distance from the helix center to linPointPos
    // in the x-y plane
    const double S = std::sqrt(S2);

    // Fill position Jacobian, i.e., matrix A from Eq. 5.36 in Ref(1)
    // First row
    positionJacobian(0, 0) = -sgnH * X / S;
    positionJacobian(0, 1) = -sgnH * Y / S;

    const double S2tanTheta = S2 * tanTheta;

    // Second row
    positionJacobian(1, 0) = rho * Y / S2tanTheta;
    positionJacobian(1, 1) = -rho * X / S2tanTheta;
    positionJacobian(1, 2) = 1.;

    // Third row
    positionJacobian(2, 0) = -Y / S2;
    positionJacobian(2, 1) = X / S2;

    // TODO: include timing in track linearization - will be added
    // in next PR
    // Sixth row
    positionJacobian(5, 3) = 1.;

    // Fill momentum Jacobian, i.e., B matrix from Eq. 5.36 in Ref(1).
    // Since we propagated to the PCA (point P in Ref(1)), the points
    // P and V coincide and we can choose deltaPhi = 0.
    // One can show that if deltaPhi = 0 --> R = 0 and Q = sgnH * S.
    // As a consequence, many terms of the B matrix from Eq. 5.36 vanish.
    double rhoOverS{sgnH * rho / S};

    // Second row
    momentumJacobian(1, 0) = (1. - rhoOverS) * rho / tanTheta;

    // Third row
    momentumJacobian(2, 0) = rhoOverS;

    // Fourth and fifth row
    momentumJacobian(3, 1) = 1.;
    momentumJacobian(4, 2) = 1.;

    // TODO: Derivatives of time --> Next PR
  }

  // const term in Talyor expansion from Eq. 5.38 in Ref(1)
  BoundVector constTerm =
      paramsAtPCA - positionJacobian * pca - momentumJacobian * momentumAtPCA;

  // The parameter weight
  BoundSymMatrix weightAtPCA = parCovarianceAtPCA.inverse();

  return LinearizedTrack(paramsAtPCA, parCovarianceAtPCA, weightAtPCA, linPoint,
                         positionJacobian, momentumJacobian, pca, momentumAtPCA,
                         constTerm);
}
