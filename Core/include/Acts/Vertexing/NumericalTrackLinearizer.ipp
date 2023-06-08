// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

template <typename propagator_t, typename propagator_options_t>
Acts::Result<Acts::LinearizedTrack>
Acts::NumericalTrackLinearizer<propagator_t, propagator_options_t>::
    linearizeTrack(const BoundTrackParameters& params, const Vector4& linPoint,
                   const Acts::GeometryContext& gctx,
                   const Acts::MagneticFieldContext& mctx) const {
  // Make Perigee surface at linPointPos, transverse plane of Perigee
  // corresponds the global x-y plane
  Vector3 linPointPos = VectorHelpers::position(linPoint);
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(linPointPos);

  // Create propagator options
  propagator_options_t pOptions(gctx, mctx);

  // Length scale at which we consider to be sufficiently close to the Perigee
  // surface to skip the propagation.
  pOptions.targetTolerance = m_cfg.targetTolerance;

  // Get intersection of the track with the Perigee if the particle would
  // move on a straight line.
  // This allows us to determine whether we need to propagate the track
  // forward or backward to arrive at the PCA.
  auto intersection = perigeeSurface->intersect(gctx, params.position(gctx),
                                                params.unitDirection(), false);

  // Setting the propagation direction using the intersection length from
  // above.
  // We handle zero path length as forward propagation, but we could actually
  // skip the whole propagation in this case.
  pOptions.direction =
      Direction::fromScalarZeroAsPositive(intersection.intersection.pathLength);

  // Propagate to the PCA of linPointPos
  auto result = m_cfg.propagator->propagate(params, *perigeeSurface, pOptions);

  // Extracting the Perigee representation of the track wrt linPointPos
  auto endParams = *result->endParameters;
  BoundVector perigeeParams = endParams.parameters();

  // Covariance and weight matrix at the PCA to "linPoint"
  BoundSymMatrix parCovarianceAtPCA = endParams.covariance().value();
  BoundSymMatrix weightAtPCA = parCovarianceAtPCA.inverse();

  // Vector containing the track parameters at the PCA
  // Note that we parametrize the track using the following parameters:
  // (x, y, z, t, phi, theta, q/p),
  // where
  // -) (x, y, z, t) is the global 4D position of the PCA
  // -) phi and theta are the global angles of the momentum at the PCA
  // -) q/p is the charge divided by the total momentum at the PCA
  const unsigned int nParams = 7;
  Acts::ActsVector<nParams> paramVec;

  // 4D PCA and the momentum of the track at the PCA
  // These quantities will be used in the computation of the constant term in
  // the Taylor expansion
  Vector4 pca;
  Vector3 momentumAtPCA;

  // Fill "paramVec", "pca", and "momentumAtPCA"
  {
    Vector3 globalCoords = endParams.position(gctx);
    ActsScalar globalTime = endParams.time();
    ActsScalar phi = perigeeParams(BoundIndices::eBoundPhi);
    ActsScalar theta = perigeeParams(BoundIndices::eBoundTheta);
    ActsScalar qOvP = perigeeParams(BoundIndices::eBoundQOverP);

    paramVec << globalCoords, globalTime, phi, theta, qOvP;
    pca << globalCoords, globalTime;
    momentumAtPCA << phi, theta, qOvP;
  }

  // Complete Jacobian (consists of positionJacobian and momentumJacobian)
  ActsMatrix<eBoundSize, nParams> completeJacobian;
  completeJacobian.setZero();

  // Perigee parameters wrt linPoint after wiggling
  BoundVector newPerigeeParams;

  // Check if wiggled angle theta are within definition range [0, pi]
  assert(paramVec(5) + m_cfg.delta < M_PI &&
         "Wiggled theta outside range, choose a smaller wiggle (i.e., delta)!"
         "You might need to decrease targetTolerance as well.");
  // Wiggling each of the parameters at the PCA and computing the Perigee
  // parametrization of the resulting new track. This allows us to approximate
  // the numerical derivatives.
  for (unsigned int i = 0; i < nParams; i++) {
    // Wiggle
    paramVec(i) += m_cfg.delta;

    // Create curvilinear track object from our parameters. This is needed for
    // the propagation. Note that we work without covariance since we don't need
    // it to compute the derivative.
    Vector3 wiggledDir =
        makeDirectionUnitFromPhiTheta(paramVec(4), paramVec(5));
    CurvilinearTrackParameters wiggledCurvilinearParams(
        paramVec.head(4), wiggledDir, paramVec(6));

    // Obtain propagation direction
    intersection =
        perigeeSurface->intersect(gctx, paramVec.head(3), wiggledDir, false);
    pOptions.direction = Direction::fromScalarZeroAsPositive(
        intersection.intersection.pathLength);

    // Unwiggle
    paramVec(i) -= m_cfg.delta;

    // Propagate to the new PCA and extract Perigee parameters
    auto newResult = m_cfg.propagator->propagate(wiggledCurvilinearParams,
                                                 *perigeeSurface, pOptions);
    auto newEndParams = (*newResult->endParameters);
    newPerigeeParams = newEndParams.parameters();

    // Computing the numerical derivatives and filling the Jacobian
    // d_0 and z_0
    completeJacobian.block<2, 1>(0, i) =
        (newPerigeeParams.head(2) - perigeeParams.head(2)) / m_cfg.delta;
    // We need to account for the periodicity of phi (see documentiation of
    // difference_periodic)
    completeJacobian(2, i) =
        Acts::detail::difference_periodic(newPerigeeParams(2), perigeeParams(2),
                                          2 * M_PI) /
        m_cfg.delta;
    // theta, q/p, and t
    completeJacobian.block<3, 1>(3, i) =
        (newPerigeeParams.tail(3) - perigeeParams.tail(3)) / m_cfg.delta;
  }

  // Extracting positionJacobian and momentumJacobian from the complete Jacobian
  ActsMatrix<eBoundSize, 4> positionJacobian =
      completeJacobian.block<eBoundSize, 4>(0, 0);
  ActsMatrix<eBoundSize, 3> momentumJacobian =
      completeJacobian.block<eBoundSize, 3>(0, 4);

  // Constant term of Taylor expansion (Eq. 5.38 in Ref. (1))
  BoundVector constTerm =
      perigeeParams - positionJacobian * pca - momentumJacobian * momentumAtPCA;

  return LinearizedTrack(perigeeParams, parCovarianceAtPCA, weightAtPCA,
                         linPoint, positionJacobian, momentumJacobian, pca,
                         momentumAtPCA, constTerm);
}
