// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterTypes.hpp"

#include <cmath>
#include <type_traits>

// The user can override the (track) parameter ordering and underlying scalar
// type. If the variable is defined, it must point to a header file that
// contains the same enum and type definitions for bound and free track
// parameters as well as space points as given below.
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
  // The following aliases without prefix exist for historical reasons
  // Generic spatial coordinates on the local surface
  eLOC_0 = eBoundLoc0,
  eLOC_1 = eBoundLoc1,
  // Spatial coordinates on a disk in polar coordinates
  eLOC_R = eLOC_0,
  eLOC_PHI = eLOC_1,
  // Spatial coordinates on a disk in Cartesian coordinates
  eLOC_X = eLOC_0,
  eLOC_Y = eLOC_1,
  // Spatial coordinates on a cylinder
  eLOC_RPHI = eLOC_0,
  eLOC_Z = eLOC_1,
  // Closest approach coordinates on a virtual perigee surface
  eLOC_D0 = eLOC_0,
  eLOC_Z0 = eLOC_1,
  // Direction angles
  ePHI = eBoundPhi,
  eTHETA = eBoundTheta,
  eQOP = eBoundQOverP,
  eT = eBoundTime,
};

/// Underlying fundamental scalar type for bound track parameters.
using BoundScalar = double;

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

/// Underlying fundamental scalar type for free track parameters.
using FreeScalar = double;

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

// Ensure bound track parameter components/ indices are consistently defined.
static_assert(eLOC_0 != eLOC_1, "Local parameters must be differents");
static_assert(eLOC_R == eLOC_0 or eLOC_R == eLOC_1,
              "Local radius must be a local parameter");
static_assert(eLOC_PHI == eLOC_0 or eLOC_PHI == eLOC_1,
              "Local phi must be a local parameter");
static_assert(eLOC_RPHI == eLOC_0 or eLOC_RPHI == eLOC_1,
              "Local r*phi must be a local parameter");
static_assert(eLOC_Z == eLOC_0 or eLOC_Z == eLOC_1,
              "Local z must be a local parameter");
static_assert(eLOC_X == eLOC_0 or eLOC_X == eLOC_1,
              "Local x must be a local parameter");
static_assert(eLOC_Y == eLOC_0 or eLOC_Y == eLOC_1,
              "Local y must be a local parameter");
static_assert(eLOC_D0 == eLOC_0 or eLOC_D0 == eLOC_1,
              "D0 must be a local parameter");
static_assert(eLOC_Z0 == eLOC_0 or eLOC_Z0 == eLOC_1,
              "Z0 must be a local parameter");

// Ensure free track parameter components/ indices are consistently defined.
static_assert(eFreePos1 == eFreePos0 + 1u, "Position must be continous");
static_assert(eFreePos2 == eFreePos0 + 2u, "Position must be continous");
static_assert(eFreeDir1 == eFreeDir0 + 1u, "Direction must be continous");
static_assert(eFreeDir2 == eFreeDir0 + 2u, "Direction must be continous");

namespace detail {

template <BoundIndices>
struct BoundParameterTraits;
template <>
struct BoundParameterTraits<BoundIndices::eBoundLoc0> {
  using type = local_parameter;
};
template <>
struct BoundParameterTraits<BoundIndices::eBoundLoc1> {
  using type = local_parameter;
};
template <>
struct BoundParameterTraits<BoundIndices::eBoundPhi> {
  static constexpr double pMin() { return -M_PI; }
  static constexpr double pMax() { return M_PI; }
  using type = cyclic_parameter<double, pMin, pMax>;
};
template <>
struct BoundParameterTraits<BoundIndices::eBoundTheta> {
  static constexpr double pMin() { return 0; }
  static constexpr double pMax() { return M_PI; }
  using type = bound_parameter<double, pMin, pMax>;
};
template <>
struct BoundParameterTraits<BoundIndices::eBoundQOverP> {
  using type = unbound_parameter;
};
template <>
struct BoundParameterTraits<BoundIndices::eBoundTime> {
  using type = unbound_parameter;
};

template <FreeIndices>
struct FreeParameterTraits;
template <>
struct FreeParameterTraits<FreeIndices::eFreePos0> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreePos1> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreePos2> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreeTime> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreeDir0> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreeDir1> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreeDir2> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreeQOverP> {
  using type = unbound_parameter;
};

template <typename indices_t>
struct ParametersSize;
template <>
struct ParametersSize<BoundIndices> {
  static constexpr unsigned int size =
      static_cast<unsigned int>(BoundIndices::eBoundSize);
};
template <>
struct ParametersSize<FreeIndices> {
  static constexpr unsigned int size =
      static_cast<unsigned int>(FreeIndices::eFreeSize);
};

}  // namespace detail

/// Single bound track parameter type for value constrains.
///
/// The singular name is not a typo since this describes individual components.
template <BoundIndices kIndex>
using BoundParameterType = typename detail::BoundParameterTraits<kIndex>::type;

/// Single free track parameter type for value constrains.
///
/// The singular name is not a typo since this describes individual components.
template <FreeIndices kIndex>
using FreeParameterType = typename detail::FreeParameterTraits<kIndex>::type;

/// Access the (Bound/Free)ParameterType through common struct
///
/// @tparam T Parameter indices enum
/// @tparam I Index of @p T
template <typename T, T I>
struct ParameterTypeFor {};

/// Access for @c BoundIndices
template <BoundIndices I>
struct ParameterTypeFor<BoundIndices, I> {
  using type = BoundParameterType<I>;
};

/// Access for @c FreeIndices
template <FreeIndices I>
struct ParameterTypeFor<FreeIndices, I> {
  using type = FreeParameterType<I>;
};

// The following matrix and vector types are automatically derived from the
// indices enums and scalar typedefs.

// Matrix and vector types related to bound track parameters.

using BoundVector = ActsVector<BoundScalar, eBoundSize>;
using BoundRowVector = ActsRowVector<BoundScalar, eBoundSize>;
using BoundMatrix = ActsMatrix<BoundScalar, eBoundSize, eBoundSize>;
using BoundSymMatrix = ActsSymMatrix<BoundScalar, eBoundSize>;

using LocalCartesianToBoundLocalMatrix = ActsMatrix<BoundScalar, 2, 3>;

// Matrix and vector types related to free track parameters.

using FreeVector = ActsVector<FreeScalar, eFreeSize>;
using FreeRowVector = ActsRowVector<FreeScalar, eFreeSize>;
using FreeMatrix = ActsMatrix<FreeScalar, eFreeSize, eFreeSize>;
using FreeSymMatrix = ActsSymMatrix<FreeScalar, eFreeSize>;

// Mapping to bound track parameters.
//
// Assumes that matrices represent maps from another space into the space of
// bound track parameters. Thus, the bound scalar type is sufficient
// to retain accuracy.

using FreeToBoundMatrix = ActsMatrix<BoundScalar, eBoundSize, eFreeSize>;

// Mapping to free track parameters.
//
// Assumes that matrices represent maps from another space into the space of
// free track parameters. Thus, the free scalar type is sufficient
// to retain accuracy.

using BoundToFreeMatrix = ActsMatrix<FreeScalar, eFreeSize, eBoundSize>;

// For backward compatibility. New code must use the more explicit
// `Bound{Indices,Scalar,Traits}...` types.
using ParDef = BoundIndices;
using ParID_t = BoundIndices;
using ParValue_t = BoundScalar;

}  // namespace Acts
