// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PerigeeSurface.hpp"

template <typename propagator_t, typename propagator_options_t>
Acts::Result<Acts::LinearizedTrack>
Acts::NumericalTrackLinearizer<propagator_t, propagator_options_t>::
    linearizeTrack(const BoundTrackParameters& params, const Vector4& linPoint,
                   const Acts::GeometryContext& gctx,
                   const Acts::MagneticFieldContext& mctx) const {
  // Make Perigee surface at linPointPos, transverse plane of Perigee
  // corresponds the global x-y plane
  Vector3 linPointPos{VectorHelpers::position(linPoint)};
  std::shared_ptr<PerigeeSurface> perigeeSurface{
      Surface::makeShared<PerigeeSurface>(linPointPos)};

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
  BoundSymMatrix parCovarianceAtPCA{endParams.covariance().value()};
  BoundSymMatrix weightAtPCA{parCovarianceAtPCA.inverse()};

  // Vector containing curvilinear track parameters at the PCA
  const unsigned int nCurvilinearParams{7};
  Acts::ActsVector<nCurvilinearParams> curvilinearParamVec;

  // 4D PCA and the momentum of the track at the PCA
  // These quantities will be used in the computation of the constant term in
  // the Taylor expansion
  Vector4 pca;
  Vector3 momentumAtPCA;

  // Fill "curvilinearParamVec", "pca", and "momentumAtPCA"
  {
    Vector3 globalCoords{endParams.position(gctx)};
    double globalTime{endParams.time()};
    double phi{perigeeParams(BoundIndices::eBoundPhi)};
    double theta{perigeeParams(BoundIndices::eBoundTheta)};
    double qOvP{perigeeParams(BoundIndices::eBoundQOverP)};

    curvilinearParamVec << globalCoords, globalTime, phi, theta, qOvP;
    pca << globalCoords, globalTime;
    momentumAtPCA << phi, theta, qOvP;
  }

  // Setting size of the perturbation delta for calculation of numerical
  // derivatives (i.e., f'(x) ~ (f(x+delta) - f(x)) / delta)
  double delta{1e-8};

  // Complete Jacobian (consists of positionJacobian and momentumJacobian)
  ActsMatrix<eBoundSize, nCurvilinearParams> completeJacobian;
  completeJacobian.setZero();

  // Curvilinear parameters at the PCA to linPoint after wiggling
  Vector4 newPos;
  Vector3 newDir;
  ActsScalar newP;
  ActsScalar newQ;

  // Perigee parameters wrt linPoint after wiggling
  BoundVector newPerigeeParams;

  // Wiggling each of the curvilinear parameters at the PCA and computing the
  // Perigee parametrization of the resulting new track. This allows us to
  // approximate the numerical derivatives.
  for (unsigned int i = 0; i < nCurvilinearParams; i++) {
    // Wiggle
    curvilinearParamVec(i) += delta;
    newPos = curvilinearParamVec.head(4);
    newDir = Vector3(
        std::sin(curvilinearParamVec(5)) * std::cos(curvilinearParamVec(4)),
        std::sin(curvilinearParamVec(5)) * std::sin(curvilinearParamVec(4)),
        std::cos(curvilinearParamVec(5)));
    newP = std::abs(1.0 / curvilinearParamVec(6));
    newQ = (curvilinearParamVec(6) > 0 ? 1.0 : -1.0);
    curvilinearParamVec(i) -= delta;

    // Curvilinear parameters object needed for the propagation
    CurvilinearTrackParameters newCurvilinearParams(newPos, newDir, newP, newQ);

    // Obtain propagation direction
    intersection =
        perigeeSurface->intersect(gctx, newPos.head(3), newDir, false);
    pOptions.direction = Direction::fromScalarZeroAsPositive(
        intersection.intersection.pathLength);

    // Propagate to the new PCA and extract Perigee parameters
    auto newResult = m_cfg.propagator->propagate(newCurvilinearParams,
                                                 *perigeeSurface, pOptions);
    auto newEndParams = (*newResult->endParameters);
    newPerigeeParams = newEndParams.parameters();

    // Computing the numerical derivatives and filling the Jacobian
    completeJacobian.array().col(i) =
        (newPerigeeParams - perigeeParams) / delta;
  }

  // Extracting positionJacobian and momentumJacobian from the complete Jacobian
  ActsMatrix<eBoundSize, 4> positionJacobian{
      completeJacobian.block<eBoundSize, 4>(0, 0)};
  ActsMatrix<eBoundSize, 3> momentumJacobian{
      completeJacobian.block<eBoundSize, 3>(0, 4)};

  // Constant term of Taylor expansion (Eq. 5.38 in Ref. (1))
  BoundVector constTerm{perigeeParams - positionJacobian * pca -
                        momentumJacobian * momentumAtPCA};

  return LinearizedTrack(perigeeParams, parCovarianceAtPCA, weightAtPCA,
                         linPoint, positionJacobian, momentumJacobian, pca,
                         momentumAtPCA, constTerm);
}
