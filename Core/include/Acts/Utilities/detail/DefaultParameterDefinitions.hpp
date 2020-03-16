// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

#include "Acts/Utilities/ParameterTypes.hpp"

namespace Acts {

/// Components of a bound track parameters vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
enum BoundParametersIndices : unsigned int {
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
  eBoundParametersSize,
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
  BoundParsDim = eBoundParametersSize,
};

/// Underlying fundamental scalar type for bound track parameters.
using BoundParametersScalar = double;

/// Components of a free track parameters vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
enum FreeParametersIndices : unsigned int {
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
  // See BoundParametersIndices for further information
  eFreeQOverP = 7u,
  // Last uninitialized value contains the total number of components
  eFreeParametersSize,
  // For backward compatibility
  FreeParsDim = eFreeParametersSize,
};

/// Underlying fundamental scalar type for free track parameters.
using FreeParametersScalar = double;

/// Components of a space point vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
///
/// Within the same context either the position-like or the momentum-like
/// indices must be used exclusively.
enum SpacePointIndices : unsigned int {
  // For position four-vectors
  // The spatial position components must be stored as one continous block.
  eSpacePos0 = 0u,
  eSpacePos1 = eSpacePos0 + 1u,
  eSpacePos2 = eSpacePos0 + 2u,
  eSpaceTime = 3u,
  // Last uninitialized value contains the total number of components
  eSpacePointSize,
  // Aliases for  momentum four-vectors to allow clearer code
  eSpaceMom0 = eSpacePos0,
  eSpaceMom1 = eSpacePos1,
  eSpaceMom2 = eSpacePos2,
  eSpaceEnergy = eSpaceTime,
  // for backward compatibility
  SpacePointDim = eSpacePointSize,
};

/// Underlying fundamental scalar type for space points.
using SpacePointScalar = double;

using ParDef = BoundParametersIndices;
using ParID_t = BoundParametersIndices;
using ParValue_t = BoundParametersScalar;

template <ParID_t>
struct par_type;

template <ParID_t par>
using par_type_t = typename par_type<par>::type;

template <>
struct par_type<ParDef::eLOC_0> {
  using type = local_parameter;
};

template <>
struct par_type<ParDef::eLOC_1> {
  using type = local_parameter;
};

template <>
struct par_type<ParDef::ePHI> {
  static constexpr double pMin() { return -M_PI; }
  static constexpr double pMax() { return M_PI; }
  using type = cyclic_parameter<double, pMin, pMax>;
};

template <>
struct par_type<ParDef::eTHETA> {
  static constexpr double pMin() { return 0; }
  static constexpr double pMax() { return M_PI; }
  using type = bound_parameter<double, pMin, pMax>;
};

template <>
struct par_type<ParDef::eQOP> {
  using type = unbound_parameter;
};

template <>
struct par_type<ParDef::eT> {
  using type = unbound_parameter;
};

}  // namespace Acts
