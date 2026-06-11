// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/Helpers.hpp"

#include <iostream>

#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h>

namespace ActsPlugins::ActsToMille {

Acts::DynamicMatrix regulariseCovariance(const Acts::DynamicMatrix& inputCov,
                                         double conditionCutOff,
                                         double removeHugeLeading,
                                         double stabilisationDiag) {
  Acts::DynamicMatrix out =
      Acts::DynamicMatrix::Zero(inputCov.rows(), inputCov.cols());

  /// optionally add a tiny diagonal matrix for additional stabilisation
  Acts::DynamicMatrix regularisation =
      Acts::DynamicMatrix::Zero(inputCov.rows(), inputCov.cols());
  if (stabilisationDiag > 0) {
    regularisation = stabilisationDiag * Acts::DynamicMatrix::Identity(
                                             inputCov.rows(), inputCov.cols());
  }
  auto eigensolver = Eigen::SelfAdjointEigenSolver<Acts::DynamicMatrix>(
      inputCov + regularisation);

  if (eigensolver.info() != Eigen::Success) {
    std::cout << " FAILED to find eigenvec" << std::endl;
    return out;
  }
  auto eigenVals = eigensolver.eigenvalues();
  auto eigenVecs = eigensolver.eigenvectors();

  std::size_t maxIndex = eigenVals.size() - 1;
  double lambdaMax = eigenVals(eigenVals.size() - 1);
  // check for a huge leading eigenvalue - this happens when the time coordinate
  // is unconstrained
  if (removeHugeLeading > 0 &&
      lambdaMax > removeHugeLeading * eigenVals(maxIndex - 1)) {
    --maxIndex;
    lambdaMax = eigenVals(maxIndex);
  }
  double lambdaMin = eigenVals(0);
  if (conditionCutOff > 0) {
    lambdaMin = conditionCutOff * lambdaMax;
  }
  // clamp the EV to the permitted interval
  for (Eigen::Index i = 0; i < eigenVals.size(); ++i) {
    out(i, i) = std::clamp(eigenVals(i), lambdaMin, lambdaMax);
  }
  // then return the regularised matrix as V D V^T
  out = out * eigenVecs.transpose();
  out = eigenVecs * out;
  return out;
}

Acts::DynamicMatrix getInverseComplement(
    const Acts::DynamicMatrix& target,
    const Acts::DynamicMatrix& existing_sol) {
  Acts::DynamicMatrix Rhs =
      Acts::DynamicMatrix::Identity(target.rows(), target.cols()) -
      target * existing_sol;
  // call decomposition from Eigen (prefer over llt for semi-def. matrices)
  auto LDLT = target.ldlt();
  Acts::DynamicMatrix solution = LDLT.solve(Rhs);
  // finally, symmetrise to correct for floating point effects
  solution = 0.5 * (solution + solution.transpose());
  return solution;
}

}  // end namespace ActsPlugins::ActsToMille
