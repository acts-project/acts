// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <type_traits>

#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

/// Components of alignment parameters vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
enum AlignmentParametersIndices : unsigned int {
  // Center of geometry object in global 3D cartesian coordinates
  eAlignmentCenter0 = 0u,
  eAlignmentCenter1 = eAlignmentCenter0 + 1u,
  eAlignmentCenter2 = eAlignmentCenter0 + 2u,
  // Rotation angle around global x/y/z axis of geometry object
  eAlignmentRotation0 = 3u,
  eAlignmentRotation1 = eAlignmentRotation0 + 1u,
  eAlignmentRotation2 = eAlignmentRotation0 + 2u,
  // Last uninitialized value contains the total number of components
  eAlignmentParametersSize,
};

/// Underlying fundamental scalar type for alignment parameters.
using AlignmentParametersScalar = double;

// Matrix and vector types related to alignment parameters.
using AlignmentVector =
    ActsVector<AlignmentParametersScalar, eAlignmentParametersSize>;
using AlignmentRowVector =
    ActsRowVector<AlignmentParametersScalar, eAlignmentParametersSize>;
using AlingmentMatrix =
    ActsMatrix<AlignmentParametersScalar, eAlignmentParametersSize,
               eAlignmentParametersSize>;
using AlignmentToLocalCartesianMatrix =
    ActsMatrix<AlignmentParametersScalar, 3, eAlignmentParametersSize>;
using AlignmentToBoundMatrix =
    ActsMatrix<BoundParametersScalar, eBoundParametersSize,
               eAlignmentParametersSize>;
}  // namespace Acts
