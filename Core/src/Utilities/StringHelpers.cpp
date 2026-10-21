// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/StringHelpers.hpp"

#include <sstream>

#include <Eigen/Eigenvalues>

namespace Acts {

std::string printEigenDecomposition(
    const Eigen::Ref<const Eigen::MatrixXd>& mat) {
  // The matrix is assumed to be symmetric (e.g. a Hessian), so the
  // self-adjoint solver yields real eigen values and vectors and is
  // considerably cheaper to compile than the general Eigen::EigenSolver.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenDecomp{mat};

  std::stringstream sstr{};
  sstr << "eigen values: " << eigenDecomp.eigenvalues().transpose();
  sstr << ", eigen vectors:\n" << eigenDecomp.eigenvectors();
  return sstr.str();
}

}  // namespace Acts
