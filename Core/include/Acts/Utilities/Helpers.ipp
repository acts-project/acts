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
  static_assert(!A.IsRowMajor, "not col major");
  A(0, 0) += 1;
  A(1, 1) += 1;
  A(2, 2) += 1;
  A(3, 3) += 1;
  A(4, 4) += 1;
  A(5, 5) += 1;
  A(6, 6) += 1;
  A(7, 7) += 1;

  return A;
}

}  // namespace MatrixHelpers
}  // namespace Acts
