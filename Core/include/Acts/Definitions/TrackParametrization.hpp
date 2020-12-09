// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <type_traits>

// The user can override the (track) parameter ordering and underlying
// Scalar type. If the variable is defined, it must point to a header
// file that contains the same enum and type definitions for bound and free
// track parameters as well as space points as given below.
#ifdef ACTS_PARAMETER_DEFINITIONS_HEADER
#include ACTS_PARAMETER_DEFINITIONS_HEADER
#else

namespace Acts {

// Note:
// The named indices are use to access raw data vectors and matrices at the
// lowest level. Since the interpretation of some of the components, e.g. local
// position and the inverse-momentum-like component, depend on additional
// information the names have some ambiguity. This can only be resolved at a
// higher logical level and no attempt is made to resolve it here.

/// Components of a bound track parameters vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
enum BoundIndices : unsigned int {
  // Local position on the reference surface.
  // This is intentionally named different from the position components in
  // the other data vectors, to clarify that this is defined on a surface
  // while the others are defined in free space.
  eBoundLoc0 = 0,
  eBoundLoc1 = 1,
  // Direction angles
  eBoundPhi = 2,
  eBoundTheta = 3,
  // Global inverse-momentum-like parameter, i.e. q/p or 1/p
  // The naming is inconsistent for the case of neutral track parameters where
  // the value is interpreted as 1/p not as q/p. This is intentional to avoid
  // having multiple aliases for the same element and for lack of an acceptable
  // common name.
  eBoundQOverP = 4,
  eBoundTime = 5,
  // Last uninitialized value contains the total number of components
  eBoundSize,
};

/// Underlying fundamental Scalar type for bound track parameters.
using BoundScalar = ActsScalar;

/// Components of a free track parameters vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
enum FreeIndices : unsigned int {
  // Spatial position
  // The spatial position components must be stored as one continous block.
  eFreePos0 = 0u,
  eFreePos1 = eFreePos0 + 1u,
  eFreePos2 = eFreePos0 + 2u,
  // Time
  eFreeTime = 3u,
  // (Unit) direction
  // The direction components must be stored as one continous block.
  eFreeDir0 = 4u,
  eFreeDir1 = eFreeDir0 + 1u,
  eFreeDir2 = eFreeDir0 + 2u,
  // Global inverse-momentum-like parameter, i.e. q/p or 1/p
  // See BoundIndices for further information
  eFreeQOverP = 7u,
  // Last uninitialized value contains the total number of components
  eFreeSize,
};

/// Underlying fundamental Scalar type for free track parameters.
using FreeScalar = ActsScalar;

}  // namespace Acts
#endif

namespace Acts {

// Ensure bound track parameters definition is valid.
static_assert(std::is_enum_v<BoundIndices>,
              "'BoundIndices' must be an enum type");
static_assert(std::is_convertible_v<BoundIndices, size_t>,
              "'BoundIndices' must be convertible to size_t");
static_assert(2 <= BoundIndices::eBoundSize,
              "Bound track parameters must have at least two components");
static_assert(std::is_floating_point_v<BoundScalar>,
              "'BoundScalar' must be a floating point type");

// Ensure free track parameters definition is valid.
static_assert(std::is_enum_v<FreeIndices>,
              "'FreeIndices' must be an enum type");
static_assert(std::is_convertible_v<FreeIndices, size_t>,
              "'FreeIndices' must be convertible to size_t");
static_assert(6 <= FreeIndices::eFreeSize,
              "Free track parameters must have at least six components");
static_assert(std::is_floating_point_v<FreeScalar>,
              "'FreeScalar' must be a floating point type");

// Ensure bound track parameter indices are consistently defined.
static_assert(eBoundLoc0 != eBoundLoc1, "Local parameters must be differents");

// Ensure free track parameter indices are consistently defined.
static_assert(eFreePos1 == eFreePos0 + 1u, "Position must be continous");
static_assert(eFreePos2 == eFreePos0 + 2u, "Position must be continous");
static_assert(eFreeDir1 == eFreeDir0 + 1u, "Direction must be continous");
static_assert(eFreeDir2 == eFreeDir0 + 2u, "Direction must be continous");

// The following matrix and vector types are automatically derived from the
// indices enums and Scalar typedefs.

// Matrix and vector types related to bound track parameters.

using BoundVector = Eigen::Matrix<BoundScalar, eBoundSize, 1>;
using BoundMatrix = Eigen::Matrix<BoundScalar, eBoundSize, eBoundSize>;
using BoundSymMatrix = Eigen::Matrix<BoundScalar, eBoundSize, eBoundSize>;

// Matrix and vector types related to free track parameters.

using FreeVector = Eigen::Matrix<FreeScalar, eFreeSize, 1>;
using FreeMatrix = Eigen::Matrix<FreeScalar, eFreeSize, eFreeSize>;
using FreeSymMatrix = Eigen::Matrix<FreeScalar, eFreeSize, eFreeSize>;

// Mapping from bound track parameters.
//
// Assumes that matrices represent maps from the bound parameters space into
// another space. To retain the accuracy of the input space, the mappings must
// use the same scalar as the input space.

using BoundToFreeMatrix = Eigen::Matrix<BoundScalar, eFreeSize, eBoundSize>;

// Mapping from free track parameters.
//
// Assumes that matrices represent maps from the free parameters space into
// another space. To retain the accuracy of the input space, the mappings must
// use the same scalar as the input space.

using FreeToBoundMatrix = Eigen::Matrix<FreeScalar, eBoundSize, eFreeSize>;
using FreeToPathMatrix = Eigen::Matrix<FreeScalar, 1, eFreeSize>;

// Other mappings.

// TODO this is not track parameter specific and might fit better to the general
//   algebra definitions as it connects three- and two-dimensional spaces.
//   maybe just replace w/ ActsMatrix<2,3>?
using PositionToBoundLocalMatrix = Eigen::Matrix<BoundScalar, 2, 3>;

}  // namespace Acts
