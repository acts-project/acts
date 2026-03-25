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

  SourceLinks = 1 << 0,           ///< Source link information
  X = 1 << 1,                     ///< X coordinate
  Y = 1 << 2,                     ///< Y coordinate
  Z = 1 << 3,                     ///< Z coordinate
  R = 1 << 4,                     ///< Radial coordinate
  Phi = 1 << 5,                   ///< Azimuthal angle
  Time = 1 << 6,                  ///< Time information
  VarianceZ = 1 << 7,             ///< Variance in Z direction
  VarianceR = 1 << 8,             ///< Variance in radial direction
  TopStripVector = 1 << 9,        ///< Vector for the top strip
  BottomStripVector = 1 << 10,    ///< Vector for the bottom strip
  StripCenterDistance = 1 << 11,  ///< Distance to the strip center
  TopStripCenter = 1 << 12,       ///< Center of the top strip
  CopyFromIndex = 1 << 13,        ///< Copy from index

  // packed columns for performance reasons
  PackedXY = 1 << 14,          ///< X and Y coordinates
  PackedZR = 1 << 15,          ///< Z and R coordinates
  PackedXYZ = 1 << 16,         ///< X, Y, and Z coordinates
  PackedXYZR = 1 << 17,        ///< X, Y, Z, and R coordinates
  PackedVarianceZR = 1 << 18,  ///< Variance in Z and R directions

  /// All strip-related columns
  Strip =
      TopStripVector | BottomStripVector | StripCenterDistance | TopStripCenter,

  /// All columns
  All = SourceLinks | X | Y | Z | R | Phi | Time | VarianceZ | VarianceR |
        TopStripVector | BottomStripVector | StripCenterDistance |
        TopStripCenter | CopyFromIndex | PackedXY | PackedZR | PackedXYZ |
        PackedXYZR | PackedVarianceZR,

  XY [[deprecated("Use PackedXY instead")]] = 1 << 14,
  ZR [[deprecated("Use PackedZR instead")]] = 1 << 15,
  XYZ [[deprecated("Use PackedXYZ instead")]] = 1 << 16,
  XYZR [[deprecated("Use PackedXYZR instead")]] = 1 << 17,
  VarianceZR [[deprecated("Use PackedVarianceZR instead")]] = 1 << 18,
};

/// Enable bitwise operators for SpacePointColumns enum
ACTS_DEFINE_ENUM_BITWISE_OPERATORS(SpacePointColumns);

}  // namespace Acts
