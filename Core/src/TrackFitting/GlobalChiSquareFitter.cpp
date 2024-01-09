// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GlobalChiSquareFitter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"

namespace Acts::Experimental {

BoundVector calculateDeltaParams(bool zeroField, const BoundMatrix& aMatrix,
                                 const BoundVector& bVector) {
  BoundVector deltaParams = BoundVector::Zero();
  if (zeroField) {
    constexpr std::size_t reducedMatrixSize = 4;
    deltaParams.topLeftCorner<reducedMatrixSize, 1>() =
        aMatrix.topLeftCorner<reducedMatrixSize, reducedMatrixSize>()
            .colPivHouseholderQr()
            .solve(bVector.topLeftCorner<reducedMatrixSize, 1>());
  } else {
    constexpr std::size_t reducedMatrixSize = 5;
    deltaParams.topLeftCorner<reducedMatrixSize, 1>() =
        aMatrix.topLeftCorner<reducedMatrixSize, reducedMatrixSize>()
            .colPivHouseholderQr()
            .solve(bVector.topLeftCorner<reducedMatrixSize, 1>());
  }

  return deltaParams;
}
}  // namespace Acts::Experimental
