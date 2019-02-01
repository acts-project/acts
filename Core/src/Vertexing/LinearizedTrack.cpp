// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/LinearizedTrack.hpp"

Acts::LinearizedTrack::LinearizedTrack(
    const ActsVectorD<5>&    paramsAtPCA,
    const ActsSymMatrixD<5>& parCovarianceAtPCA,
    const Vector3D&          linPoint,
    const ActsMatrixD<5, 3>& positionJacobian,
    const ActsMatrixD<5, 3>& momentumJacobian,
    const Vector3D&       positionAtPCA,
    const Vector3D&       momentumAtPCA,
    const ActsVectorD<5>& constTerm)
  : parametersAtPCA(paramsAtPCA)
  , covarianceAtPCA(parCovarianceAtPCA)
  , linearizationPoint(linPoint)
  , positionJacobian(positionJacobian)
  , momentumJacobian(momentumJacobian)
  , positionAtPCA(positionAtPCA)
  , momentumAtPCA(momentumAtPCA)
  , constantTerm(constTerm)
{
}