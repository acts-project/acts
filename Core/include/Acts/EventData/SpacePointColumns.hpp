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

  // strip-specific columns
  InnerStripCenter = 1 << 10,  ///< Center of the inner strip
  InnerStripVector = 1 << 11,  ///< Vector of the inner strip
  OuterStripCenter = 1 << 12,  ///< Center of the outer strip
  OuterStripVector = 1 << 13,  ///< Vector of the outer strip

  // derived strip columns
  /// Inner to outer strip center vector
  StripCenterDistance = 1 << 14,
  /// Cross product of strip center distance and outer strip vector
  StripCenterDistanceCrossOuterStripVector = 1 << 15,
  /// Cross product of strip center distance and inner strip vector
  StripCenterDistanceCrossInnerStripVector = 1 << 16,
  /// Cross product of inner and outer strip vectors
  InnerStripVectorCrossOuterStripVector = 1 << 17,

  // packed columns for performance reasons
  PackedXY = 1 << 18,          ///< X and Y coordinates
  PackedZR = 1 << 19,          ///< Z and R coordinates
  PackedXYZ = 1 << 20,         ///< X, Y, and Z coordinates
  PackedXYZR = 1 << 21,        ///< X, Y, Z, and R coordinates
  PackedVarianceZR = 1 << 22,  ///< Variance in Z and R directions

  /// Relevant strip columns for seeding purposes
  StripRelevant = InnerStripVector | InnerStripCenter |
                  StripCenterDistanceCrossInnerStripVector |
                  StripCenterDistanceCrossOuterStripVector |
                  InnerStripVectorCrossOuterStripVector,

  /// All strip columns
  StripAll = InnerStripVector | InnerStripCenter | OuterStripVector |
             OuterStripCenter | StripCenterDistance |
             StripCenterDistanceCrossInnerStripVector |
             StripCenterDistanceCrossOuterStripVector |
             InnerStripVectorCrossOuterStripVector,

  /// All columns
  All = SourceLinks | X | Y | Z | R | Phi | Time | VarianceZ | VarianceR |
        CopiedFromIndex | InnerStripCenter | InnerStripVector |
        OuterStripCenter | OuterStripVector | StripCenterDistance |
        StripCenterDistanceCrossInnerStripVector |
        StripCenterDistanceCrossOuterStripVector |
        InnerStripVectorCrossOuterStripVector | PackedXY | PackedZR |
        PackedXYZ | PackedXYZR | PackedVarianceZR,
};

/// Enable bitwise operators for SpacePointColumns enum
ACTS_DEFINE_ENUM_BITWISE_OPERATORS(SpacePointColumns);

}  // namespace Acts
