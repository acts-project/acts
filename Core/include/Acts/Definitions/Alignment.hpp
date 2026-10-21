// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"

namespace Acts {

/// Components of alignment parameters vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
enum AlignmentIndices : unsigned int {
  // Center of geometry object in global 3D cartesian coordinates
  eAlignmentCenter0 = 0u,
  eAlignmentCenter1 = eAlignmentCenter0 + 1u,
  eAlignmentCenter2 = eAlignmentCenter0 + 2u,
  // Rotation angle around local x/y/z axis of geometry object
  eAlignmentRotation0 = 3u,
  eAlignmentRotation1 = eAlignmentRotation0 + 1u,
  eAlignmentRotation2 = eAlignmentRotation0 + 2u,
  // Last uninitialized value contains the total number of components
  eAlignmentSize,
};

// Matrix and vector types related to alignment parameters.
/// @brief Vector type for alignment parameters
using AlignmentVector = Vector<eAlignmentSize>;
/// @brief Row vector type for alignment parameters
using AlignmentRowVector = Matrix<1, eAlignmentSize>;
/// @brief Square matrix type for alignment parameters
using AlignmentMatrix = Matrix<eAlignmentSize, eAlignmentSize>;
/// @brief Matrix type for transforming alignment parameters to position
using AlignmentToPositionMatrix = Matrix<3, eAlignmentSize>;
/// @brief Matrix type for transforming alignment parameters to bound parameters
using AlignmentToBoundMatrix = Matrix<eBoundSize, eAlignmentSize>;
/// @brief Matrix type for transforming alignment parameters to path length
using AlignmentToPathMatrix = Matrix<1, eAlignmentSize>;

}  // namespace Acts
