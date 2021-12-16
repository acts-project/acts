// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"

namespace Acts {
namespace MatrixHelpers {

#include "Acts/Utilities/HelpersMatrixGenerated.ipp"

inline ActsMatrix<8, 8> plusIdentity(ActsMatrix<8, 8> A) {
  double* pA = A.data();

  pA[0] += 1;
  pA[9] += 1;
  pA[18] += 1;
  pA[27] += 1;
  pA[36] += 1;
  pA[45] += 1;
  pA[54] += 1;
  pA[63] += 1;

  return A;
}

}  // namespace MatrixHelpers
}  // namespace Acts
