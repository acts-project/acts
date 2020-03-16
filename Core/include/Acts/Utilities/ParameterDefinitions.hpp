// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <type_traits>

#include "Acts/Utilities/Definitions.hpp"

#ifdef ACTS_PARAMETER_DEFINITIONS_PLUGIN
#include ACTS_PARAMETER_DEFINITIONS_PLUGIN
#endif

namespace Acts {

// Ensure bound track parameters definition is valid.
static_assert(std::is_enum_v<BoundParametersIndices>,
              "'BoundParametersIndices' must be an enum type");
static_assert(std::is_convertible_v<BoundParametersIndices, size_t>,
              "'BoundParametersIndices' must be convertible to size_t");
static_assert(2 <= BoundParametersIndices::eBoundParametersSize,
              "Bound track parameters must have at least two components");
static_assert(std::is_floating_point_v<BoundParametersScalar>,
              "'BoundParametersScalar' must be a floating point type");

// Ensure free track parameters definition is valid.
static_assert(std::is_enum_v<FreeParametersIndices>,
              "'FreeParametersIndices' must be an enum type");
static_assert(std::is_convertible_v<FreeParametersIndices, size_t>,
              "'FreeParametersIndices' must be convertible to size_t");
static_assert(6 <= FreeParametersIndices::eFreeParametersSize,
              "Free track parameters must have at least six components");
static_assert(std::is_floating_point_v<FreeParametersScalar>,
              "'FreeParametersScalar' must be a floating point type");

// Ensure space point definition is valid.
static_assert(std::is_enum_v<SpacePointIndices>,
              "'SpacePointIndices' is not an enum type");
static_assert(std::is_convertible_v<SpacePointIndices, size_t>,
              "'SpacePointIndices' is not convertible to size_t");
static_assert(3 <= SpacePointIndices::eSpacePointSize,
              "Space points must have at least three components");
static_assert(std::is_floating_point_v<SpacePointScalar>,
              "'SpacePointScalar' must be a floating point type");

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

// Ensure space point components/ indices are consistently defined.
static_assert(eSpacePos1 == eSpacePos0 + 1u, "Position must be continous");
static_assert(eSpacePos2 == eSpacePos0 + 2u, "Position must be continous");
static_assert(eSpacePos0 == eSpaceMom0, "Inconsisten position and momentum");
static_assert(eSpacePos1 == eSpaceMom1, "Inconsisten position and momentum");
static_assert(eSpacePos2 == eSpaceMom2, "Inconsisten position and momentum");
static_assert(eSpaceTime == eSpaceEnergy, "Inconsistent time and energy");

namespace detail {
template <BoundParametersIndices>
struct BoundParameterTraits;
template <>
struct BoundParameterTraits<BoundParametersIndices::eBoundLoc0> {
  using type = local_parameter;
};
template <>
struct BoundParameterTraits<BoundParametersIndices::eBoundLoc1> {
  using type = local_parameter;
};
template <>
struct BoundParameterTraits<BoundParametersIndices::eBoundPhi> {
  static constexpr double pMin() { return -M_PI; }
  static constexpr double pMax() { return M_PI; }
  using type = cyclic_parameter<double, pMin, pMax>;
};
template <>
struct BoundParameterTraits<BoundParametersIndices::eBoundTheta> {
  static constexpr double pMin() { return 0; }
  static constexpr double pMax() { return M_PI; }
  using type = bound_parameter<double, pMin, pMax>;
};
template <>
struct BoundParameterTraits<BoundParametersIndices::eBoundQOverP> {
  using type = unbound_parameter;
};
template <>
struct BoundParameterTraits<BoundParametersIndices::eBoundTime> {
  using type = unbound_parameter;
};
}  // namespace detail

/// Single bound track parameter type for value constrains.
///
/// The singular name is not a typo since this describes individual components.
template <BoundParametersIndices kIndex>
using BoundParameterType = typename detail::BoundParameterTraits<kIndex>::type;

// The following matrix and vector types are automatically derived from the
// indices enums and scalar typedefs.

// Matrix and vector types related to bound track parameters.

using BoundVector = ActsVector<BoundParametersScalar, eBoundParametersSize>;
using BoundRowVector =
    ActsRowVector<BoundParametersScalar, eBoundParametersSize>;
using BoundMatrix = ActsMatrix<BoundParametersScalar, eBoundParametersSize,
                               eBoundParametersSize>;
using BoundSymMatrix =
    ActsSymMatrix<BoundParametersScalar, eBoundParametersSize>;

// Matrix and vector types related to free track parameters.

using FreeVector = ActsVector<FreeParametersScalar, eFreeParametersSize>;
using FreeRowVector = ActsRowVector<FreeParametersScalar, eFreeParametersSize>;
using FreeMatrix =
    ActsMatrix<FreeParametersScalar, eFreeParametersSize, eFreeParametersSize>;
using FreeSymMatrix = ActsSymMatrix<FreeParametersScalar, eFreeParametersSize>;

// Matrix and vector types related to space points.

using SpacePointVector = ActsVector<SpacePointScalar, eSpacePointSize>;
using SpacePointRowVector = ActsRowVector<SpacePointScalar, eSpacePointSize>;
using SpacePointSymMatrix =
    ActsMatrix<SpacePointScalar, eSpacePointSize, eSpacePointSize>;
using SpacePointSymMatrix = ActsSymMatrix<SpacePointScalar, eSpacePointSize>;

// Mapping to bound track parameters.
//
// Assumes that matrices represent maps from another space into the space of
// bound track parameters. Thus, the bound parameters scalar type is sufficient
// to retain accuracy.

using FreeToBoundMatrix = ActsMatrix<BoundParametersScalar,
                                     eBoundParametersSize, eFreeParametersSize>;
using SpacePointToBoundMatrix =
    ActsMatrix<BoundParametersScalar, eBoundParametersSize, eSpacePointSize>;

// Mapping to free track parameters.
//
// Assumes that matrices represent maps from another space into the space of
// free track parameters. Thus, the free parameters scalar type is sufficient
// to retain accuracy.

using BoundToFreeMatrix =
    ActsMatrix<FreeParametersScalar, eFreeParametersSize, eBoundParametersSize>;
using SpacePointToFreeMatrix =
    ActsMatrix<FreeParametersScalar, eFreeParametersSize, eSpacePointSize>;

// Mapping to space points.
//
// Assumes that matrices represent maps from another space into the space point
// space. Thus, the space point scalar type is sufficient to retain accuracy.

using BoundToSpacePointMatrix =
    ActsMatrix<SpacePointScalar, eSpacePointSize, eBoundParametersSize>;
using FreeToSpacePointMatrix =
    ActsMatrix<SpacePointScalar, eSpacePointSize, eFreeParametersSize>;

// For backward compatibility. New code must use the more explicit
// `BoundParameters{Indices,Scalar,Traits}...` types.
using ParDef = BoundParametersIndices;
using ParID_t = BoundParametersIndices;
using ParValue_t = BoundParametersScalar;
template <BoundParametersIndices kIndex>
using par_type = detail::BoundParameterTraits<kIndex>;
template <BoundParametersIndices kIndex>
using par_type_t = typename detail::BoundParameterTraits<kIndex>::type;

// For backward compatibility. New code must use the
// `e{BoundParameters,FreeParameters,SpacePoint}Size` enum values.
inline constexpr unsigned int BoundParsDim = eBoundParametersSize;
inline constexpr unsigned int FreeParsDim = eFreeParametersSize;
inline constexpr unsigned int SpacePointDim = eSpacePointSize;

}  // namespace Acts
