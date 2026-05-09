// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/EnumBitwiseOperators.hpp"

#include <cstdint>

namespace Acts {

/// Enumeration of available columns for space point data storage
enum class SpacePointColumns : std::uint32_t {
  None = 0,  ///< No columns

  SourceLinks = 1 << 0,      ///< Source link information
  CopiedFromIndex = 1 << 1,  ///< Copied from index
  X = 1 << 2,                ///< X coordinate
  Y = 1 << 3,                ///< Y coordinate
  Z = 1 << 4,                ///< Z coordinate
  R = 1 << 5,                ///< Radial coordinate
  Phi = 1 << 6,              ///< Azimuthal angle
  Time = 1 << 7,             ///< Time information
  VarianceZ = 1 << 8,        ///< Variance in Z direction
  VarianceR = 1 << 9,        ///< Variance in radial direction

  /// Calibration details for strip space points
  StripCalibrationDetails = 1 << 10,

  // packed columns for performance reasons
  PackedXY = 1 << 11,          ///< X and Y coordinates
  PackedZR = 1 << 12,          ///< Z and R coordinates
  PackedXYZ = 1 << 13,         ///< X, Y, and Z coordinates
  PackedXYZR = 1 << 14,        ///< X, Y, Z, and R coordinates
  PackedVarianceZR = 1 << 15,  ///< Variance in Z and R directions

  /// All columns
  All = SourceLinks | X | Y | Z | R | Phi | Time | VarianceZ | VarianceR |
        CopiedFromIndex | StripCalibrationDetails | PackedXY | PackedZR |
        PackedXYZ | PackedXYZR | PackedVarianceZR,
};

/// Enable bitwise operators for SpacePointColumns enum
ACTS_DEFINE_ENUM_BITWISE_OPERATORS(SpacePointColumns);

}  // namespace Acts
